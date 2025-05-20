#include <fstream>
#include <iostream>
#include <string>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <experiments.hpp>

#include <chrono>
#include <mockturtle/algorithms/reconv_cut.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/utils/window_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/window_view.hpp>

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

template<typename Ntk>
kitty::dynamic_truth_table compute_care_set( Ntk ntk, uint32_t index, const std::vector<mockturtle::node<Ntk>>& leaves )
{
  static constexpr uint32_t window_size = 14;
  static constexpr uint32_t max_window_size = 14;

  reconvergence_driven_cut_parameters rps;
  rps.max_leaves = window_size; // ps.window_size;
  reconvergence_driven_cut_statistics rst;
  detail::reconvergence_driven_cut_impl<Ntk, false, false> reconv_cuts( ntk, rps, rst );

  mockturtle::color_view<Ntk> color_ntk{ ntk };

  std::vector<mockturtle::node<Ntk>> roots = {  index };
  auto const extended_leaves = reconv_cuts.run( roots ).first;

  auto depth_ntk = mockturtle::depth_view(mockturtle::fanout_view(color_ntk));
  create_window_impl windowing(depth_ntk);
  const auto res = windowing.run(roots[0], window_size, 5);

  if ( res != std::nullopt )
  {
    int zz = 0;
  }

  std::vector<mockturtle::node<Ntk>> gates{ collect_nodes( color_ntk, extended_leaves, roots ) };
  window_view window_ntk{ color_ntk, extended_leaves, roots, gates };

  if ( leaves.size() == extended_leaves.size() )
  {
    std::cout << "Happens\n";
    kitty::dynamic_truth_table care = kitty::dynamic_truth_table( leaves.size() );
    return ~care;
  }
  std::cout << "Leaves size: " << leaves.size() << ", Extended Leaves size: " << extended_leaves.size() << std::endl;

  // dont cares computation
  kitty::dynamic_truth_table care = kitty::dynamic_truth_table( leaves.size() );

  bool containment = true;
  for ( auto const& l : leaves )
  {
    if ( color_ntk.color( ntk.index_to_node( l ) ) != color_ntk.current_color() )
    {
      containment = false;
      break;
    }
  }

  if ( containment )
  {
    default_simulator<kitty::static_truth_table<max_window_size>> sim;
    const auto tts = simulate_nodes<kitty::static_truth_table<max_window_size>>( window_ntk, sim );
    // compute care set
    for ( auto i = 0u; i < ( 1u << window_ntk.num_pis() ); ++i )
    {
      uint32_t entry{ 0u };
      auto j = 0u;
      for ( auto const& l : leaves )
      {
        const auto lr = tts[l];
        entry |= kitty::get_bit( tts[l], i ) << j;
        ++j;
      }
      kitty::set_bit( care, entry );
    }
  }
  else
  {
    care = ~care;
  }
  return care;
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
    /*if ( benchmark != "sin" )
    {
      continue;
    }*/
    /*if ( benchmark != "mem_ctrl" )
    {
      continue;
    }*/
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    std::string bench_path = "/home/benjamin/Documents/Project_dump/ACD/ACD_benchmarks/opt_benchmarks_temp/";
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
    uint32_t num_dcs = 0;
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
        std::vector<mockturtle::node<aig_network>> leaves;
        for (const auto& l : *cut)
        {
          leaves.push_back(l);
        }
        const auto care = compute_care_set( aig, index, leaves );

        uint32_t const num_blocks = ( cut_size > 6 ) ? ( 1u << ( cut_size - 6 ) ) : 1;
        for ( uint32_t i = 0; i < num_blocks; ++i )
        {
          const auto c = care._bits[i];
          num_dcs += static_cast<uint32_t>(64 - std::__popcount(c)); // count unset bits
        }

        /*for ( uint32_t i = 0; i < num_blocks; ++i )
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
        }*/
        ++num_cuts;
      }
    } );
    std::cout << "Num Cuts: " << num_cuts << std::endl;
    std::cout << "Num DCs: " << num_dcs << std::endl;
    // break;
    exp( benchmark, aig.num_pis(), aig.size(), num_cuts );

    break;
  }
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Execution time: " << elapsed_seconds.count() << "s\n";

  exp.save();
  exp.table();

  return 0;
}