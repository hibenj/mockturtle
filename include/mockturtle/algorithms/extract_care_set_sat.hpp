/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
 \file extract_care_sets_sat.hpp
 \brief Generate miter from two networks

 \author Benjamin Hien
*/

#pragma once

#include <cassert>
#include <cstdint>
#include <vector>
#include <bill/sat/solver/abc.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <mockturtle/algorithms/cnf.hpp>

namespace mockturtle
{

namespace detail
{

template<typename Ntk, typename Cut>
class extract_care_set_sat_impl
{
public:
  extract_care_set_sat_impl( Ntk const& ntk, Cut const* const& cut )
      : ntk_( ntk ),
        cut_( cut ),
        global_lits_( node_literals( ntk ) ),
        solver_( pabc::sat_solver_new() )
  {
    assert( cut_ != nullptr );

    int output_lit = generate_cnf( ntk, [&]( auto const& clause ) {
      auto lits = (int*)( const_cast<uint32_t*>( clause.data() ) );
      pabc::sat_solver_addclause( solver_, lits, lits + clause.size() );
    } )[0];

    output_lit_ = output_lit;
  }

  ~extract_care_set_sat_impl()
  {
    pabc::sat_solver_delete( solver_ );
  }

  kitty::dynamic_truth_table run()
  {
    const uint32_t num_leaves = cut_->size();
    kitty::dynamic_truth_table care_tt( num_leaves );

    std::vector<uint32_t> cut_leaves;
    for ( auto leaf : *cut_ )
      cut_leaves.insert( cut_leaves.begin(), leaf );

    for ( int polarity = 0; polarity <= 1; ++polarity )
    {
      int assumption = ( polarity == 0 ) ? output_lit_ : ( output_lit_ ^ 1 );

      while ( true )
      {
        int status = pabc::sat_solver_solve( solver_, &assumption, &assumption + 1, 0, 0, 0, 0 );
        if ( status == -1 )
          break;

        std::vector<uint32_t> cut_blocking_clause;
        // std::vector<uint32_t> pi_blocking_clause;
        uint64_t index = 0;

        // 1. Block cut leaf assignment
        for ( size_t i = 0; i < cut_leaves.size(); ++i )
        {
          const auto node = cut_leaves[i];
          const auto lit = global_lits_[node];
          bool value = pabc::sat_solver_var_value( solver_, node );
          cut_blocking_clause.push_back( value ? lit ^ 1 : lit );

          // Also build the index for the truth table
          index |= ( uint64_t( value ) << ( cut_leaves.size() - 1 - i ) );
        }
        // Block full PI assignment too
        /*ntk_.foreach_pi([&](const auto& pi) {
          const auto lit = global_lits_[pi];
          bool value = pabc::sat_solver_var_value(solver_, pi);
          pi_blocking_clause.push_back(value ? lit ^ 1 : lit);  // block this PI assignment
        });*/

        // Set the corresponding bit in the truth table
        care_tt._bits[index >> 6] |= uint64_t( 1 ) << ( index & 0x3F );

        // Add both blocking clauses
        auto cut_lits = (int*)( cut_blocking_clause.data() );
        pabc::sat_solver_addclause( solver_, cut_lits, cut_lits + cut_blocking_clause.size() );

        /*auto pi_lits = (int*)(pi_blocking_clause.data());
        pabc::sat_solver_addclause(solver_, pi_lits, pi_lits + pi_blocking_clause.size());*/
      }
    }

    return care_tt;
  }

private:
  Ntk const& ntk_;
  Cut const* const& cut_;
  pabc::sat_solver* solver_;
  node_map<uint32_t, Ntk> global_lits_;
  int output_lit_;
};

// Public interface
template<typename Ntk, typename Cut>
kitty::dynamic_truth_table extract_care_set_sat( Ntk const& ntk, Cut const* const& cut )
{
  detail::extract_care_set_sat_impl<Ntk, Cut> impl( ntk, cut );
  return impl.run();
}

template<typename Ntk>
class extract_care_set_sat_from_leaves_impl
{
public:
  extract_care_set_sat_from_leaves_impl( Ntk const& ntk, std::vector<uint32_t> const& leaves )
      : ntk_( ntk ),
        leaves_( leaves ),
        global_lits_( node_literals( ntk ) ),
        solver_( pabc::sat_solver_new() )
  {
    int output_lit = generate_cnf( ntk, [&]( auto const& clause ) {
      auto lits = (int*)( const_cast<uint32_t*>( clause.data() ) );
      pabc::sat_solver_addclause( solver_, lits, lits + clause.size() );
    } )[0];

    output_lit_ = output_lit;
  }

  ~extract_care_set_sat_from_leaves_impl()
  {
    pabc::sat_solver_delete( solver_ );
  }

  kitty::dynamic_truth_table run()
  {
    kitty::dynamic_truth_table care_tt( leaves_.size() );

    for ( int polarity = 0; polarity <= 1; ++polarity )
    {
      int assumption = ( polarity == 0 ) ? output_lit_ : ( output_lit_ ^ 1 );

      while ( true )
      {
        int status = pabc::sat_solver_solve( solver_, &assumption, &assumption + 1, 0, 0, 0, 0 );
        if ( status == -1 )
          break;

        std::vector<uint32_t> cut_blocking_clause;
        uint64_t index = 0;

        for ( size_t i = 0; i < leaves_.size(); ++i )
        {
          const auto node = leaves_[i];
          const auto lit = global_lits_[node];
          bool value = pabc::sat_solver_var_value( solver_, node );
          cut_blocking_clause.push_back( value ? lit ^ 1 : lit );

          index |= ( uint64_t( value ) << ( leaves_.size() - 1 - i ) );
        }

        care_tt._bits[index >> 6] |= uint64_t( 1 ) << ( index & 0x3F );

        auto cut_lits = (int*)( cut_blocking_clause.data() );
        pabc::sat_solver_addclause( solver_, cut_lits, cut_lits + cut_blocking_clause.size() );
      }
    }

    return care_tt;
  }

private:
  Ntk const& ntk_;
  std::vector<uint32_t> leaves_;
  pabc::sat_solver* solver_;
  node_map<uint32_t, Ntk> global_lits_;
  int output_lit_;
};

// Public interface
template<typename Ntk>
kitty::dynamic_truth_table extract_care_set_sat( Ntk const& ntk, std::vector<uint32_t> const& leaves )
{
  detail::extract_care_set_sat_from_leaves_impl<Ntk> impl( ntk, leaves );
  return impl.run();
}

} // namespace detail

} // namespace mockturtle
