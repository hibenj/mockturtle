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
#include <mockturtle/algorithms/dont_cares.hpp>
#include <mockturtle/algorithms/extract_care_set_sat.hpp>
#include <mockturtle/algorithms/reconv_cut.hpp>

using namespace experiments;
using namespace mockturtle;

int main()
{
  experiment<std::string, uint32_t, uint32_t, double, uint32_t, uint32_t, double, float, float, bool, bool> exp(
      "mapper", "benchmark", "size", "size_mig", "area_after", "depth", "depth_mig", "delay_after", "runtime1", "runtime2", "equivalent1", "equivalent2" );

  for ( auto& benchmark : epfl_benchmarks() )
  {
    if ( benchmark != "ctrl" )
    {
      continue;
    }
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    std::string bench_path = "/home/benjamin/Desktop/Notes/ACD_benchmarks/opt_benchmarks_temp/";
    std::string bench = bench_path + benchmark + ".aig";
    if ( lorina::read_aiger( bench, aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    static constexpr unsigned cut_size = 7u;
    cut_enumeration_params ps_c;
    ps_c.cut_size = cut_size;
    ps_c.cut_limit = 40u;
    cut_enumeration_stats st_c;
    const auto aig_topo = mockturtle::topo_view( aig );
    const auto cuts = fast_cut_enumeration<decltype( aig_topo ), cut_size, true, cut_enumeration_params>( aig_topo, ps_c, &st_c );

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
        const auto care_set = extract_care_set_sat( aig, cut );
        std::cout << "TT:\n";
        kitty::print_binary( tt );
        std::cout << "\n";
        std::cout << "CS:\n";
        kitty::print_binary( care_set );
        std::cout << "\n";
        std::cout << kitty::count_ones( care_set ) << "\n";
        ++num_cuts;
      }
    } );
    std::cout << "Num Cuts: " << num_cuts << std::endl;
    break;
  }

  exp.save();
  exp.table();

  return 0;
}