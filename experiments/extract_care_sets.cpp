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

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>

#include <experiments.hpp>
#include <mockturtle/algorithms/reconv_cut.hpp>

using namespace experiments;
using namespace mockturtle;

void write_word_to_file( uint64_t word )
{
  std::string filename = "/home/benjamin/Desktop/Notes/ACD_benchmarks/truth_table.data";

  std::ofstream file( filename, std::ios::binary | std::ios::app ); // Open in binary append mode
  if ( !file )
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  // Write the raw 8-byte (64-bit) word directly in binary format
  file.write( reinterpret_cast<const char*>( &word ), sizeof( uint64_t ) );

  file.close();
}

template<typename Ntk, typename Cut, unsigned NInputs>
void extract_care_sets( Ntk ntk, const Cut& cuts )
{
  uint64_t num_cuts = 0;
  static constexpr uint32_t window_size = 16;
  static constexpr uint32_t max_window_size = 16;

  reconvergence_driven_cut_parameters rps;
  rps.max_leaves = window_size; // ps.window_size;
  reconvergence_driven_cut_statistics rst;
  detail::reconvergence_driven_cut_impl<Ntk, false, false> reconv_cuts( ntk, rps, rst );

  color_view<Ntk> color_ntk{ ntk };

  // match gates
  ntk.foreach_gate( [&]( auto const& n ) {
    const auto index = ntk.node_to_index( n );

    std::vector<node<Ntk>> roots = { n };
    auto const extended_leaves = reconv_cuts.run( roots ).first;

    std::vector<node<Ntk>> gates{ collect_nodes( color_ntk, extended_leaves, roots ) };
    window_view window_ntk{ color_ntk, extended_leaves, roots, gates };

    default_simulator<kitty::static_truth_table<max_window_size>> sim;
    const auto tts = simulate_nodes<kitty::static_truth_table<max_window_size>>( window_ntk, sim );

    for ( auto& cut : cuts.cuts( index ) )
    {
      // ignore unit cut
      if ( cut->size() == 1 && *cut->begin() == index )
      {
        continue;
      }

      if ( cut->size() != NInputs )
      {
        // Ignore cuts which are not maximum size
        continue;
      }

      ++num_cuts;

      // match the cut using canonization and get the gates
      const auto tt = cuts.truth_table( *cut );
      const auto fe = kitty::shrink_to<NInputs>( tt );

      // dont cares computation
      kitty::static_truth_table<NInputs> care;

      bool containment = true;
      bool filter = false;
      for ( auto const& l : *cut )
      {
        if ( color_ntk.color( ntk.index_to_node( l ) ) != color_ntk.current_color() )
        {
          containment = false;
          break;
        }
      }

      if ( containment )
      {
        // compute care set
        for ( auto i = 0u; i < ( 1u << window_ntk.num_pis() ); ++i )
        {
          uint32_t entry{ 0u };
          auto j = 0u;
          for ( auto const& l : *cut )
          {
            entry |= kitty::get_bit( tts[l], i ) << j;
            ++j;
          }
          kitty::set_bit( care, entry );
        }
      }
      else
      {
        // completely specified
        care = ~care;
      }

      uint32_t const num_blocks = ( NInputs > 6 ) ? ( 1u << ( NInputs - 6 ) ) : 1;

      if constexpr ( num_blocks == 1 )
      {
        if ( care._bits != ( 1u << ( 1u << cut->size() ) ) - 1 )
        {
          // this does not happen for var > 6
          write_word_to_file( tt._bits );
        }
      }
      else
      {
        for ( uint32_t i = 0; i < num_blocks; ++i )
        {
          // write tt to file
          const auto cof = tt._bits[i];
          write_word_to_file( cof );
        }
        for ( uint32_t i = 0; i < num_blocks; ++i )
        {
          // write cs to file
          const auto ccs = care._bits[i];
          write_word_to_file( ccs );
        }
      }
    }
  } );
  std::cout << "Num Cuts: " << num_cuts << std::endl;
}

int main()
{
  experiment<std::string, uint32_t, uint32_t, double, uint32_t, uint32_t, double, float, float, bool, bool> exp(
      "mapper", "benchmark", "size", "size_mig", "area_after", "depth", "depth_mig", "delay_after", "runtime1", "runtime2", "equivalent1", "equivalent2" );

  fmt::print( "[i] processing technology library\n" );

  for ( auto& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    mig_network mig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    aig_network aig;
    std::string bench_path = "/home/benjamin/Desktop/Notes/ACD_benchmarks/opt_benchmarks_temp/";
    std::string bench = bench_path + benchmark + ".aig";
    if ( lorina::read_aiger( bench, aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    static constexpr unsigned cut_size = 10u;
    cut_enumeration_params ps_c;
    ps_c.cut_size = cut_size;
    ps_c.cut_limit = 40u;
    cut_enumeration_stats st_c;
    const auto aig_topo = mockturtle::topo_view( aig );
    const auto cuts = fast_cut_enumeration<decltype( aig_topo ), cut_size, true, cut_enumeration_params>( aig_topo, ps_c, &st_c );

    extract_care_sets<aig_network, decltype( cuts ), cut_size>( aig, cuts );
  }

  exp.save();
  exp.table();

  return 0;
}