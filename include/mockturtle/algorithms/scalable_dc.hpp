/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
 \file window_utils.hpp
 \brief Utilities to collect small-scale sets of nodes

 \author Benjamin Hien
*/

#pragma once

#include <mockturtle/utils/window_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/window_view.hpp>
#include "simulation.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <queue>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace mockturtle
{

template<typename Ntk>
class create_dc_windows_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  struct window
  {
    std::vector<node> inputs;
    std::vector<node> nodes;
    std::vector<signal> outputs;
  };

protected:
  /* constant node used to denote invalid window elements */
  static constexpr node INVALID_NODE{ 0 };

  static constexpr uint32_t NUM_ITERATIONS{ 5 };

public:
  explicit create_dc_windows_impl( Ntk const& ntk )
      : ntk( ntk ),
        path( ntk.size() ),
        refs( ntk.size() )
  {
  }

  void resize( uint32_t size )
  {
    path.resize( size );
    refs.resize( size );
  }

  auto begin() const { return windows.begin(); }
  auto end() const { return windows.end(); }

  void mark_nodes( std::vector<node> const& leaves, uint32_t num_levels )
  {
    for ( uint32_t i = 0; i < leaves.size(); ++i )
    {
      node const& leaf = leaves[i];
      uint64_t const bit = uint64_t( 1 ) << i;

      std::queue<std::pair<node, uint32_t>> queue;
      queue.emplace( leaf, 0 );

      while ( !queue.empty() )
      {
        auto [n, level] = queue.front();
        queue.pop();

        if ( level > num_levels )
          continue;

        if ( ( visited[n] & bit ) != 0 )
          continue; // already marked for this leaf

        visited[n] |= bit;

        ntk.foreach_fanin( n, [&]( signal const& fi ) {
          queue.emplace( ntk.get_node( fi ), level + 1 );
        } );
      }
    }
  }

  bool can_expand_for_free( const node& n )
  {
    if ( ntk.is_constant( n ) || ntk.is_ci( n ) )
      return false; // constants or PIs cannot be expanded

    uint32_t count_visited_fanins = 0;

    ntk.foreach_fanin( n, [&]( auto const& fi ) {
      auto const fn = ntk.get_node( fi );
      if ( ntk.eval_color( fn, [&]( auto c ) { return c == ntk.current_color(); } ) )
        ++count_visited_fanins;
    } );

    return ( count_visited_fanins + 1 >= ntk.fanin_size( n ) );
  }

  std::vector<uint64_t> get_leaves()
  {
    std::vector<uint64_t> groups;

    std::vector<uint64_t> bitsets;
    for ( auto const& [n, bits] : visited )
    {
      assert( bits != 0 );
      bitsets.push_back( bits );
    }

    std::vector<bool> merged( bitsets.size(), false );

    for ( size_t i = 0; i < bitsets.size(); ++i )
    {
      if ( merged[i] )
        continue;

      uint64_t group = bitsets[i];
      merged[i] = true;

      size_t merge_border = bitsets.size();
      bool grew = true;
      while ( grew )
      {
        grew = false;
        size_t upper_border = merge_border;
        for ( size_t j = i + 1; j < upper_border; ++j )
        {
          if ( merged[j] )
            continue;

          if ( ( group & bitsets[j] ) != 0 )
          {
            uint64_t old_group = group;
            group |= bitsets[j];
            merged[j] = true;

            if ( std::__popcount( group ) > std::__popcount( old_group ) )
            {
              grew = true;
              merge_border = j + 1;
            }
          }
        }
      }

      groups.push_back( group );
    }

    return groups;
  }

  bool expand_cut_zero_cost( std::vector<node>& inputs )
  {
    bool changed = false;
    std::vector<node> new_inputs;

    for ( auto it = inputs.begin(); it != inputs.end(); )
    {
      if ( !can_expand_for_free( *it ) )
      {
        ++it;
        continue;
      }

      ntk.foreach_fanin( *it, [&]( auto const& fi ) {
        node fn = ntk.get_node( fi );
        if ( !ntk.eval_color( fn, [&]( auto c ) { return c == ntk.current_color(); } ) )
        {
          ntk.paint( fn );
          new_inputs.push_back( fn ); // stage new nodes safely
        }
      } );

      it = inputs.erase( it ); // remove expanded node
      changed = true;
    }

    // Safely add new nodes after loop
    inputs.insert( inputs.end(), new_inputs.begin(), new_inputs.end() );

    return changed;
  }

  bool expand_cut_information( std::vector<node>& inputs )
  {
    bool changed = false;
    std::vector<node> new_inputs;

    for ( auto it = inputs.begin(); it != inputs.end(); )
    {
      const auto visited_pop = __builtin_popcountll( visited[*it] );
      bool should_expand = false;

      std::unordered_set<node> seen;
      std::queue<node> queue;
      queue.push( *it );

      while ( !queue.empty() )
      {
        node n = queue.front();
        queue.pop();

        if ( visited.find( n ) == visited.end() || visited[n] == 0 )
          continue;

        if ( !seen.insert( n ).second )
          continue;

        if ( __builtin_popcountll( visited[n] ) > visited_pop )
        {
          should_expand = true;
          break;
        }

        ntk.foreach_fanin( n, [&]( auto const& fi ) {
          queue.push( ntk.get_node( fi ) );
        } );
      }

      if ( !should_expand )
      {
        ++it;
        continue;
      }

      // Expand node
      ntk.foreach_fanin( *it, [&]( signal const& fi ) {
        node const& fn = ntk.get_node( fi );
        if ( ntk.visited( fn ) != ntk.trav_id() )
        {
          ntk.paint( fn );
          new_inputs.push_back( fn );
        }
      } );

      it = inputs.erase( it ); // erase the expanded node
      changed = true;
    }

    // Append staged nodes after loop finishes
    inputs.insert( inputs.end(), new_inputs.begin(), new_inputs.end() );

    return changed;
  }

  void build_window( std::vector<node> const& leaves, uint64_t leave_set )
  {
    std::vector<node> inputs;
    std::vector<signal> outputs;

    ntk.new_color(); // Reset color before marking

    //std::cout << "Set\n";
    uint64_t bits = leave_set;
    while ( bits != 0 )
    {
      uint32_t i = __builtin_ctzll( bits ); // Index of the least-significant flipped bit
      bits &= bits - 1;                     // Clear the lowest flipped bit

      node const& leaf = leaves[i];
      inputs.push_back( leaf );
      outputs.push_back( ntk.make_signal( leaf ) );
      ntk.paint( leaf );
      //std::cout << leaf << std::endl;
    }
    const auto roots = inputs;

    bool changed = true;
    while ( changed )
    {
      changed = false;

      // Step 2: zero-cost expansion
      if ( expand_cut_zero_cost( inputs ) )
        changed = true;

      // Step 3: information-driven expansion
      if ( expand_cut_information( inputs ) )
        changed = true;

      if ( inputs.size() > 13 )
      {
        break;
      }
    }

    if ( changed != true )
    {
      const auto nodes = collect_nodes( ntk, inputs, roots );

      windows.push_back( window{ inputs, nodes, outputs } );
    }
  }

  void run( std::vector<node> const& leaves, uint32_t num_levels )
  {
    windows.clear(); // clear previous state

    mark_nodes( leaves, num_levels ); // traverse from the leaves through the TFI and mark as visited

    const auto window_leaves = get_leaves(); // leaves that appear at common nodes are merged into one window

    for ( const auto& leave_set : window_leaves )
    {
      build_window( leaves, leave_set ); // build a window for each independent leave_set
    }
  }

protected:
  Ntk const& ntk;
  std::vector<node> path;
  std::vector<uint32_t> refs;
  std::vector<std::vector<node>> levels;

  std::vector<window> windows;
  std::unordered_map<node, uint64_t> visited;
};

template<typename Ntk>
kitty::dynamic_truth_table scalable_dc(const Ntk& ntk, std::vector<typename Ntk::node> const& leaves, uint32_t num_levels = 12)
{
  mockturtle::color_view<Ntk> color_ntk{ntk};

  create_dc_windows_impl<decltype(color_ntk)> cut_window(color_ntk);
  cut_window.run(leaves, num_levels);

  kitty::dynamic_truth_table global_care(leaves.size());
  global_care = ~global_care;

  for (const auto& w : cut_window)
  {
    kitty::dynamic_truth_table care(w.outputs.size());
    window_view window_ntk{color_ntk, w.inputs, w.outputs, w.nodes};

    default_simulator<kitty::dynamic_truth_table> sim(w.inputs.size());
    const auto tts = simulate_nodes<kitty::dynamic_truth_table>(window_ntk, sim);

    for (auto i = 0u; i < (1u << window_ntk.num_pis()); ++i)
    {
      uint32_t entry = 0u;
      auto j = 0u;
      for (auto const& l : w.outputs)
      {
        entry |= kitty::get_bit(tts[l], i) << j;
        ++j;
      }
      kitty::set_bit(care, entry);
    }

    care = ~care;

    std::vector<kitty::dynamic_truth_table> vars;
    for (auto const& l : w.outputs)
    {
      auto it = std::find(leaves.begin(), leaves.end(), ntk.get_node(l));
      assert(it != leaves.end());
      uint32_t idx = std::distance(leaves.begin(), it);

      kitty::dynamic_truth_table var(leaves.size());
      kitty::create_nth_var(var, idx);
      vars.push_back(var);
    }

    const auto global_dc = kitty::compose_truth_table(care, vars);
    global_care &= ~global_dc;
  }

  return global_care;
}

} // namespace mockturtle