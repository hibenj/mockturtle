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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/algorithms/extract_care_set_sat.hpp>
#include <experiments.hpp>

#include <chrono>

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

void write_word_to_file_cs( uint64_t word )
{
  std::string filename = "/home/benjamin/Desktop/Notes/ACD_benchmarks/cs.data";

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

int main()
{
  experiment<std::string, uint32_t, uint32_t, uint32_t> exp( "Stats", "benchmark", "PIs", "Gates", "Cuts" );

  auto start = std::chrono::high_resolution_clock::now();
  for ( auto& benchmark : epfl_benchmarks() )
  {
    /*if ( benchmark != "div" )
    {
      continue;
    }*/
    /*if ( benchmark != "hyp" )
    {
      continue;
    }*/
    /*if ( benchmark != "log2" )
    {
      continue;
    }*/
    /*if ( benchmark != "multiplier" )
    {
      continue;
    }*/
    if ( benchmark != "sin" )
    {
      continue;
    }
    /*if ( benchmark != "mem_ctrl" )
    {
      continue;
    }*/
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    std::string bench_path = "/home/benjamin/Desktop/Notes/ACD_benchmarks/opt_benchmarks_temp/";
    std::string bench = bench_path + benchmark + ".aig";
    if ( lorina::read_aiger( bench, aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    /*if( aig.num_gates() > 16000 )
    {
      continue;
    }*/

    static constexpr unsigned cut_size = 7u;
    cut_enumeration_params ps_c;
    ps_c.cut_size = cut_size;
    ps_c.cut_limit = 40u;
    cut_enumeration_stats st_c;
    const auto aig_topo = mockturtle::topo_view( aig );
    const auto cuts = fast_cut_enumeration<decltype( aig_topo ), cut_size, true, cut_enumeration_params>( aig_topo, ps_c, &st_c );
    fmt::print( "[i] cut enumeration finished {}\n", benchmark );
    uint32_t num_cuts = 0;
    aig_topo.foreach_gate( [&]( auto const& n ) {
      const auto index = aig_topo.node_to_index( n );
      for ( auto& cut : cuts.cuts( index ) )
      {
        if ( cut->size() != cut_size )
        {
          // Ignore cuts which are not maximum size
          continue;
        }

        const auto tt = cuts.truth_table( *cut );
        const auto care = mockturtle::detail::extract_care_set_sat( aig, cut );

        uint32_t const num_blocks = ( cut_size > 6 ) ? ( 1u << ( cut_size - 6 ) ) : 1;
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
          write_word_to_file_cs( ccs );
        }
        std::cout << "Num Cuts: " << num_cuts << std::endl;
        ++num_cuts;
      }
    } );
    std::cout << "Num Cuts: " << num_cuts << std::endl;
    // break;
    exp( benchmark, aig.num_pis(), aig.size(), num_cuts );
  }
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Execution time: " << elapsed_seconds.count() << "s\n";

  exp.save();
  exp.table();

  return 0;
}