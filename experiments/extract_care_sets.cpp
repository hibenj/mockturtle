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
#include <mockturtle/algorithms/acd.hpp>
#include <mockturtle/algorithms/acd_dc.hpp>

#include <experiments.hpp>
#include <mockturtle/algorithms/dont_cares.hpp>
// #include <mockturtle/algorithms/extract_care_set_sat.hpp>
#include <mockturtle/algorithms/reconv_cut.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>

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

template<typename Ntk, typename Cut, unsigned NInputs>
void extract_care_sets( const Ntk& ntk, const Cut& cuts, std::string benchmark )
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
      continue;

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
        /*if ( benchmark == "hyp" )
        {
          continue;
        }
        if ( benchmark == "log2" )
        {
          continue;
        }
        if ( benchmark == "multiplier" )
        {
          continue;
        }
        if ( benchmark == "sqrt" )
        {
          continue;
        }
        if ( benchmark == "mem_ctrl" )
        {
          continue;
        }*/

        // const auto care_dyn = extract_care_set_sat( ntk, cut );
        // care = care_dyn;
        care = ~care;
      }
      // ++num_cuts;
      // std::cout << "Num Cuts: " << benchmark << " progress " << num_cuts << std::endl;
      /*std::cout << "TT:\n";
      kitty::print_binary(tt);
      std::cout << "\n";
      std::cout << "CS:\n";
      kitty::print_binary(care);
      std::cout << "\n";
      std::cout << kitty::count_ones(care) << "\n";*/

      // write to file
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
          write_word_to_file_cs( ccs );
        }
      }
    }
  } );
  std::cout << "Num Cuts: " << num_cuts << std::endl;
}

uint64_t binomial_coefficient( uint32_t n, uint32_t k )
{
  if ( k > n ) return 0;
  if ( k == 0 || k == n ) return 1;

  uint64_t result = 1;
  for ( uint32_t i = 1; i <= k; ++i )
  {
    result = result * ( n - i + 1 ) / i;
  }
  return result;
}

int32_t mockturtle_acd_generic( std::string const& tt_string, uint32_t delay_profile )
{
  using namespace mockturtle;

  uint32_t num_vars = std::log2( 4 * tt_string.size() );

  kitty::dynamic_truth_table tt( num_vars );
  kitty::create_from_hex_string( tt, tt_string );

  uint64_t words[1024];

  for ( uint32_t i = 0; i < tt.num_blocks(); ++i )
    words[i] = tt._bits[i];

  acd_params ps;
  ps.use_first = false;
  ps.max_multiplicity = 16;
  acd_stats st;
  acd_impl acd( num_vars, ps, &st );

  int res = acd.run( words, delay_profile );

  if ( res < 0 )
    return -1;

  return st.num_luts;
}

int32_t mockturtle_acd_dc_generic( std::string const& tt_string, std::string const& cs_string, uint32_t delay_profile )
{
  using namespace mockturtle;

  uint32_t num_vars = std::log2( 4 * tt_string.size() );
  uint32_t num_vars_cs = std::log2( 4 * cs_string.size() );

  kitty::dynamic_truth_table tt( num_vars );
  kitty::create_from_hex_string( tt, tt_string );

  kitty::dynamic_truth_table cs( num_vars );
  kitty::create_from_hex_string( cs, cs_string );

  uint64_t words[1024];
  uint64_t words_cs[1024];

  for ( uint32_t i = 0; i < tt.num_blocks(); ++i )
    words[i] = tt._bits[i];
  for ( uint32_t i = 0; i < cs.num_blocks(); ++i )
    words_cs[i] = cs._bits[i];

  acd_params ps;
  ps.use_first = false;
  ps.max_multiplicity = 16;
  acd_stats st;
  acd_dc_impl acd_dc( num_vars, ps, &st );

  int res = acd_dc.run( words, words_cs, delay_profile );

  if ( res < 0 )
    return -1;

  return st.num_luts;
}

void compute_success_rate_delay( uint32_t cut_size, uint32_t late_vars = 2 )
{
  /* read file */
  std::string in_string = "cuts_" + std::to_string( cut_size ) + ".txt";
  std::string in_cs_string = "cuts_cs_" + std::to_string( cut_size ) + ".txt";
  // std::string in_string = "unopt/var7/cuts_" + std::to_string( cut_size ) + "_16.txt";
  // std::string in_cs_string = "unopt/var7/cuts_cs_" + std::to_string( cut_size ) + "_16.txt";
  std::ifstream in( in_string );
  std::ifstream in_cs( in_cs_string );

  /* count the number of lines */
  uint32_t num_lines = 0;
  std::string line;
  while ( std::getline( in, line ) )
  {
    ++num_lines;
  }

  in.close();

  /* count the number of lines */
  uint32_t num_lines_cs = 0;
  std::string line_cs;
  while ( std::getline( in_cs, line_cs ) )
  {
    ++num_lines_cs;
  }

  in_cs.close();

  assert( num_lines == num_lines_cs );

  in.open( in_string );
  in_cs.open( in_cs_string );

  std::ofstream out( "cuts_" + std::to_string( cut_size ) + "_fail.txt" );

  using clock = typename std::chrono::steady_clock;
  using time_point = typename clock::time_point;

  time_point time_begin = clock::now();
  std::default_random_engine::result_type seed{ 1 };
  std::uniform_int_distribution<uint32_t> dist( 0u, cut_size - 1 );

  /* compute */
  uint32_t successJ = 0, successG = 0;
  uint64_t num_LutsJ = 0, num_LutsG=0;
  uint32_t visit = 0;
  std::chrono::duration<double> total_time_generic{0};
  std::chrono::duration<double> total_time_acd_dc{0};
  while ( in.good() )
  {
    std::cout << fmt::format( "[i] Progress {:8d} / {}\r", visit, num_lines );
    std::string tt;
    in >> tt;
    std::string cs;
    in_cs >> cs;

    ++visit;

    if ( tt.size() < 16 )
      continue;

    /* generate random delay profile with late variables */

    uint32_t delay = (1u << late_vars) - 1;  // first bitmask with late_vars bits set
    // delay = 0b10010010101;
    // std::cout << std::bitset<32>(delay) << std::endl;
    while (true)
    {
      uint32_t delay_profile = delay;

      /* run evaluation */
      // Time ACD_DC version
      time_point start_dc = clock::now();
      int32_t resJ = mockturtle_acd_dc_generic( tt, cs, delay_profile );
      total_time_acd_dc += std::chrono::duration_cast<std::chrono::duration<double>>( clock::now() - start_dc );

      // Time generic version
      time_point start_g = clock::now();
      int32_t resG = mockturtle_acd_generic( tt, delay_profile );
      total_time_generic += std::chrono::duration_cast<std::chrono::duration<double>>( clock::now() - start_g );

      if ( resG > 0 )
      {
        ++successG;
        num_LutsG += resG;
      }
      if ( resJ > 0 )
      {
        ++successJ;
        num_LutsJ += resJ;
      }

      // Generate next combination using Gosper’s hack
      uint32_t x = delay;
      if (x == ((1u << late_vars) - 1) << (cut_size - late_vars)) break;

      uint32_t smallest = x & -x;
      uint32_t ripple = x + smallest;
      uint32_t new_bits = ((ripple ^ x) >> 2) / smallest;
      delay = ripple | new_bits;
    }
  }

  std::cout << "\n";

  std::cout << "Num Luts: " << num_LutsG << std::endl;
  std::cout << "Success count: " << successG << std::endl;

  uint32_t repeat = binomial_coefficient( cut_size, late_vars );
  // uint32_t repeat = 1;

  /* print stats */
  std::cout << fmt::format( "[i] Run a total of {} truth tables on {} variables\n", num_lines, cut_size );
  std::cout << fmt::format( "[i] Success of Generic   = {} \t {:>5.2f}%\n", successG, ( (double)successG ) / num_lines * 100 / repeat );
  std::cout << fmt::format( "[i] Avg Area   = {:>5.2f}\n", ( (double)num_LutsG ) / successG );
  std::cout << fmt::format( "[i] Time (Generic)       = {:>5.2f} s\n", std::chrono::duration_cast<std::chrono::duration<double>>( total_time_generic ).count() / repeat );
  std::cout << fmt::format( "[i] Success of ACD_DC  = {} \t {:>5.2f}%\n", successJ, ( (double)successJ ) / num_lines * 100 / repeat );
  std::cout << fmt::format( "[i] Avg Area   = {:>5.2f}\n", ( (double)num_LutsJ ) / successJ );
  std::cout << fmt::format( "[i] Time (ACD_DC)        = {:>5.2f} s\n",  std::chrono::duration_cast<std::chrono::duration<double>>( total_time_acd_dc ).count() / repeat );
  std::cout << fmt::format( "[i] Time = {:>5.2f} s\n", std::chrono::duration_cast<std::chrono::duration<double>>( clock::now() - time_begin ).count() );

  in.close();
  in_cs.close();
  out.close();
}

void compute_functions( uint32_t cut_size )
{
  using namespace experiments;
  using namespace mockturtle;

  truth_table_cache<kitty::dynamic_truth_table> cache( 200000 );
  std::vector<kitty::dynamic_truth_table> cache_map;
  uint64_t size = 0;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /*aig_network aig;
    std::string bench_path = "/home/benjamin/Desktop/Notes/ACD_benchmarks/opt_benchmarks_temp/";
    std::string bench = bench_path + benchmark + ".aig";
    if ( lorina::read_aiger( bench, aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }*/

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = cut_size;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.area_share_rounds = 0;
    ps.recompute_cuts = true;
    ps.cut_expansion = false;
    lut_map_stats st;

    detail::lut_map_impl<aig_network, true, lut_unitary_cost> p( aig, ps, st );
    const auto klut = p.run();

    truth_table_cache<kitty::dynamic_truth_table> const& tt_cache = p.get_truth_cache();
    const auto& cs_cache = p.get_care_cache();

    std::cout << "Cache size: " << tt_cache.size() << std::endl;

    /* load content into cache */
    for ( uint32_t i = 0; i < tt_cache.size(); ++i )
    {
      kitty::dynamic_truth_table tt = tt_cache[i << 1];
      kitty::dynamic_truth_table cs = cs_cache[i];

      if ( tt.num_vars() != cut_size )
        continue;

      if ( cut_size <= 6 )
      {
        auto res = kitty::exact_npn_canonization( tt );
        cache.insert( std::get<0>( res ) );
      }
      else
      {
        auto res = kitty::sifting_npn_canonization( tt );
        auto res_cs = kitty::apply_npn_transformation( cs, 0, std::get<2>( res ) );
        std::size_t sz0 = cache.size();
        cache.insert( std::get<0>( res ) );
        if ( sz0 != cache.size() )
        {
          cache_map.push_back( res_cs );
        }
      }
    }
  }
  std::cout << "Cache size: " << cache.size() << std::endl;
  std::cout << "Cache_Map size: " << cache_map.size() << std::endl;

  assert( cache_map.size() == cache.size() );

  /* print truth tables to file */
  std::string filename = "cuts_" + std::to_string( cut_size ) + ".txt";
  std::ofstream out( filename );

  std::string filename_cs = "cuts_cs_" + std::to_string( cut_size ) + ".txt";
  std::ofstream out_cs( filename_cs );

  for ( uint32_t i = 0; i < cache.size(); ++i )
  {
    kitty::dynamic_truth_table tt = cache[i << 1];
    kitty::dynamic_truth_table cs = cache_map[i];

    /*for ( int j = 0; j < cs.num_blocks(); ++j )
    {
      write_word_to_file( tt._bits[j] );
    }

    for ( int j = 0; j < cs.num_blocks(); ++j )
    {
      write_word_to_file_cs( cs._bits[j] );
      write_word_to_file( cs._bits[j] );
    }*/

    kitty::print_hex( tt, out );
    out << "\n";
    kitty::print_hex( cs, out_cs );
    out_cs << "\n";
  }

  out.close();
  out_cs.close();
}

int main( int argc, char** argv )
{
  if ( argc != 2 )
    return -1;

  uint32_t cut_size = atoi( argv[1] );

  //compute_functions(cut_size);
  compute_success_rate_delay( cut_size, 5 );

  /*for ( auto& benchmark : epfl_benchmarks() )
  {
    *//*if ( benchmark != "sqrt" )
    {
      continue;
    }*//*
    fmt::print( "[i] processing {}\n", benchmark );

    // read benchmark
    aig_network aig;
    std::string bench_path = "/home/benjamin/Desktop/Notes/ACD_benchmarks/opt_benchmarks_temp/";
    std::string bench = bench_path + benchmark + ".aig";
    if ( lorina::read_aiger( bench, aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    // cut enumeration
    static constexpr unsigned cut_size = 8u;
    cut_enumeration_params ps_c;
    // ps_c.cut_size = cut_size;
    // ps_c.cut_limit = 40u;
    cut_enumeration_stats st_c;
    using CutData = cut_enumeration_tech_map_cut;
    const auto aig_topo = mockturtle::topo_view( aig );
    const auto cuts = fast_cut_enumeration<decltype( aig_topo ), cut_size, true, CutData>( aig_topo, ps_c, &st_c );

    //
    extract_care_sets<aig_network, decltype( cuts ), cut_size>( aig, cuts, benchmark );
  }*/

  return 0;
}