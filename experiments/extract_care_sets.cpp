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

#include <iostream>
#include <fstream>
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

template<typename STT>
void swap_inplace_local( STT& tt, uint8_t var_index1, uint8_t var_index2, uint32_t num_vars )
{
  if ( var_index1 == var_index2 )
  {
    return;
  }

  if ( var_index1 > var_index2 )
  {
    std::swap( var_index1, var_index2 );
  }

  const uint32_t num_blocks = num_vars <= 6 ? 1 : 1 << ( num_vars - 6 );

  if ( num_vars <= 6 )
  {
    const auto& pmask = kitty::detail::ppermutation_masks[var_index1][var_index2];
    const auto shift = ( 1 << var_index2 ) - ( 1 << var_index1 );
    tt._bits[0] = ( tt._bits[0] & pmask[0] ) | ( ( tt._bits[0] & pmask[1] ) << shift ) | ( ( tt._bits[0] & pmask[2] ) >> shift );
  }
  else if ( var_index2 <= 5 )
  {
    const auto& pmask = kitty::detail::ppermutation_masks[var_index1][var_index2];
    const auto shift = ( 1 << var_index2 ) - ( 1 << var_index1 );
    std::transform( std::begin( tt._bits ), std::begin( tt._bits ) + num_blocks, std::begin( tt._bits ),
                    [shift, &pmask]( uint64_t word ) {
                      return ( word & pmask[0] ) | ( ( word & pmask[1] ) << shift ) | ( ( word & pmask[2] ) >> shift );
                    } );
  }
  else if ( var_index1 <= 5 ) /* in this case, var_index2 > 5 */
  {
    const auto step = 1 << ( var_index2 - 6 );
    const auto shift = 1 << var_index1;
    auto it = std::begin( tt._bits );
    while ( it != std::begin( tt._bits ) + num_blocks )
    {
      for ( auto i = decltype( step ){ 0 }; i < step; ++i )
      {
        const auto low_to_high = ( *( it + i ) & kitty::detail::projections[var_index1] ) >> shift;
        const auto high_to_low = ( *( it + i + step ) << shift ) & kitty::detail::projections[var_index1];
        *( it + i ) = ( *( it + i ) & ~kitty::detail::projections[var_index1] ) | high_to_low;
        *( it + i + step ) = ( *( it + i + step ) & kitty::detail::projections[var_index1] ) | low_to_high;
      }
      it += 2 * step;
    }
  }
  else
  {
    const auto step1 = 1 << ( var_index1 - 6 );
    const auto step2 = 1 << ( var_index2 - 6 );
    auto it = std::begin( tt._bits );
    while ( it != std::begin( tt._bits ) + num_blocks )
    {
      for ( auto i = 0; i < step2; i += 2 * step1 )
      {
        for ( auto j = 0; j < step1; ++j )
        {
          std::swap( *( it + i + j + step1 ), *( it + i + j + step2 ) );
        }
      }
      it += 2 * step2;
    }
  }
}

std::string uint64_to_binary(uint64_t n)
{
  std::bitset<64> binary(n);
  return binary.to_string();
}

template<typename DTT, unsigned NInputs>
void delete_var_from_tt( uint32_t i, DTT best_tt, kitty::static_truth_table<NInputs> best_tt_cs )
{
  uint32_t constexpr num_blocks = ( NInputs > 6 ) ? ( 1u << ( NInputs - 6 ) ) : 1;
  // reduce the num_vars
  uint64_t num_vars = NInputs;

  // manipulate the truth table
  kitty::static_truth_table<NInputs> new_tt = best_tt;
  kitty::static_truth_table<NInputs> new_cs = best_tt_cs;

  auto new_tt_binary0 = uint64_to_binary(new_tt._bits[0]);
  auto new_tt_binary1 = uint64_to_binary(new_tt._bits[1]);
  auto new_tt_binary2 = uint64_to_binary(new_tt._bits[2]);
  auto new_tt_binary3 = uint64_to_binary(new_tt._bits[3]);
  auto new_cs_binary0 = uint64_to_binary(new_cs._bits[0]);
  auto new_cs_binary1 = uint64_to_binary(new_cs._bits[1]);
  auto new_cs_binary2 = uint64_to_binary(new_cs._bits[2]);
  auto new_cs_binary3 = uint64_to_binary(new_cs._bits[3]);

  // reposition the independent variable
  // ToDo: consider permutation
  for (uint32_t k = i; k < num_vars - 1; ++k) {
    swap_inplace_local(new_tt, k, k + 1, NInputs);
    swap_inplace_local(new_cs, k, k + 1, NInputs);
  }

  // halve the truth table
  uint32_t block_shift_max = 1u << ( num_vars - 7 );
  for (uint32_t k = 0; k < ( num_blocks >> 1u); ++k) {
    uint32_t l = k + block_shift_max;
    uint64_t& cof_new1 = new_tt._bits[k];
    uint64_t& ccs_new1 = new_cs._bits[k];
    uint64_t cof_new2 = new_tt._bits[l];
    uint64_t ccs_new2 = new_cs._bits[l];

    // compute the new cofactor
    cof_new1 = (cof_new1 & ccs_new1) | (cof_new2 & ccs_new2);

    // compute the new care_set
    ccs_new1 |= ccs_new2;
  }
  // set the rest of the truth table to 0
  for (uint32_t k = (num_blocks >> 1u); k < num_blocks; ++k) {
    new_tt._bits[k] = 0u;
  }

  new_tt_binary0 = uint64_to_binary(new_tt._bits[0]);
  new_tt_binary1 = uint64_to_binary(new_tt._bits[1]);
  new_tt_binary2 = uint64_to_binary(new_tt._bits[2]);
  new_tt_binary3 = uint64_to_binary(new_tt._bits[3]);
  new_cs_binary0 = uint64_to_binary(new_cs._bits[0]);
  new_cs_binary1 = uint64_to_binary(new_cs._bits[1]);
  new_cs_binary2 = uint64_to_binary(new_cs._bits[2]);
  new_cs_binary3 = uint64_to_binary(new_cs._bits[3]);
  num_vars -= 1;

  // try to reduce again
  // reduce_support();
}

template<typename DTT, unsigned NInputs>
int reduce_support( DTT best_tt, kitty::static_truth_table<NInputs> best_tt_cs )
{
  uint32_t constexpr num_blocks = ( NInputs > 6 ) ? ( 1u << ( NInputs - 6 ) ) : 1;
  kitty::static_truth_table<NInputs> v;

  for (uint32_t i = 0; i < NInputs; ++i) {
    // xs.emplace_back(push);
    kitty::create_nth_var( v, i );
    // const auto test = uint64_to_binary(v._bits[0]);

    if (i < 6) // compare inside blocks
    {
      for (uint32_t j = 0; j < num_blocks; ++j) {
        uint64_t cof;
        uint64_t ccs;
        uint64_t var;
        if constexpr ( num_blocks > 1 )
        {
          cof = best_tt._bits[j];
          ccs = best_tt_cs._bits[j];
          var = v._bits[j];
        }
        else
        {
          cof = best_tt._bits;
          ccs = best_tt_cs._bits;
          var = v._bits;
        }

        uint64_t cof1 = cof & var;
        uint64_t cof2 = ( cof & ~var) << ( 1u << i);
        uint64_t ccs1 = ccs & var;
        uint64_t ccs2 = ( ccs & ~var) << ( 1u << i);

        const auto cof_binary = uint64_to_binary(cof);
        const auto ccs_binary = uint64_to_binary(ccs);

        if ( ( cof1 ^ cof2 ) & ( ccs1 & ccs2 ) ) // the block essentially depends on variable i, therefore break and try other variable
        {
          break;
        }
        else if ( j == num_blocks - 1 ) // all blocks do not essentially depend on variable i
        {
          // delete_var_from_tt<decltype( best_tt ), NInputs>( i, best_tt, best_tt_cs );
          return 1;
        }
      }
    }
    else // compare blocks
    {
      assert(num_blocks < 33);
      uint64_t mask = ( 1u << ( num_blocks ) ) - 1 ;
      kitty::static_truth_table<7> pattern_creator;
      kitty::create_nth_var( pattern_creator, i - 6 );
      uint64_t pattern = pattern_creator._bits[0] & mask;
      uint32_t block_shift = 1u << ( i - 6 );

      while (pattern) {
        const auto pattern_binary = uint64_to_binary( pattern );
        uint32_t j = __builtin_ctzll(pattern);  // Get the lowest set bit position
        uint32_t k = j - block_shift;      // Compute k

        uint64_t cof1;
        uint64_t ccs1;
        uint64_t cof2;
        uint64_t ccs2;

        if constexpr ( num_blocks > 1 )
        {
          cof1 = best_tt._bits[j];
          ccs1 = best_tt_cs._bits[j];
          cof2 = best_tt._bits[k];
          ccs2 = best_tt_cs._bits[k];
        }
        else
        {
          assert(false);
        }

        if ( ( cof1 ^ cof2 ) & ( ccs1 & ccs2 ) ) // the block essentially depends on variable i, therefore break and try other variable
        {
          break;
        }
        else if ( j == num_blocks - 1 ) // all blocks do not essentially depend on variable i
        {
          // delete_var_from_tt<decltype( best_tt ), NInputs>( i, best_tt, best_tt_cs );
          return 1;
        }

        pattern &= pattern - 1;  // Clear the lowest set bit
      }
    }
  }
  return 0;
}

void write_word_to_file( uint64_t word )
{
  std::string filename = "/home/benjamin/Desktop/Notes/ACD_benchmarks/truth_table.data";

  std::ofstream file(filename, std::ios::binary | std::ios::app); // Open in binary append mode
  if (!file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  // Write the raw 8-byte (64-bit) word directly in binary format
  file.write(reinterpret_cast<const char*>(&word), sizeof(uint64_t));

  file.close();
}

template< typename Ntk, typename Cut, unsigned NInputs>
void extract_care_sets(Ntk ntk, const Cut& cuts, uint32_t& reduce)
{
  uint64_t num_cuts = 0;
  static constexpr uint32_t window_size = 12;
  static constexpr uint32_t max_window_size = 12;

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
      /*kitty::static_truth_table<NInputs> cs;
      for ( uint32_t i = 0; i < num_blocks; ++i )
      {
        cs._bits[i] = 0xFFFFFFFFFFFFFFFE;
      }*/
      // reduce += reduce_support<decltype( tt ), NInputs>( tt, care );

      if constexpr (num_blocks == 1)
      {
        if ( care._bits != ( 1u << ( 1u << cut->size() ) ) - 1 )
        {
          // write_word_to_file( tt._bits );
        }
      }
      else
      {
        for (uint32_t i = 0; i < num_blocks; ++i)
        {
          const auto cof = tt._bits[i]; // Use iterator if more than one block
          const auto ccs = care._bits[i];
          write_word_to_file( cof );
          write_word_to_file( ccs );
        }
      }

      int zzu = 0;
    }
  } );
  std::cout << "Num Cuts: " << num_cuts  << std::endl;
}

int main()
{
 experiment<std::string, uint32_t, uint32_t, double, uint32_t, uint32_t, double, float, float, bool, bool> exp(
     "mapper", "benchmark", "size", "size_mig", "area_after", "depth", "depth_mig", "delay_after", "runtime1", "runtime2", "equivalent1", "equivalent2" );

 fmt::print( "[i] processing technology library\n" );

 for ( auto& benchmark : epfl_benchmarks() )
 {
   if (benchmark != "max" )
   {
     continue;
   }
   /*else
   {
     benchmark = "priority_opt";
   }*/
   fmt::print( "[i] processing {}\n", benchmark );
   mig_network mig;
   if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) ) != lorina::return_code::success )
   {
     continue;
   }

   aig_network aig;
   if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
   {
     continue;
   }

   static constexpr unsigned cut_size = 7u;
   cut_enumeration_params ps_c;
   ps_c.cut_size = cut_size;
   ps_c.cut_limit = 40u;
   cut_enumeration_stats st_c;
   const auto aig_topo = mockturtle::topo_view(aig);
   const auto cuts = fast_cut_enumeration<decltype(aig_topo), cut_size, true, cut_enumeration_params>( aig_topo, ps_c, &st_c );
   // const auto cuts = cut_enumeration<aig_network, true>( aig_topo, ps_c, &st_c );

   uint32_t reduce = 0;
   extract_care_sets<aig_network, decltype( cuts ), cut_size>(aig, cuts, reduce);

   // std::cout << "num_cuts: " << num_cuts << std::endl;
   // std::cout << "reduce: " << reduce << std::endl;

   // exp( benchmark, size_before, res1.num_gates(), st2.area, depth_before, depth_mig, st2.delay, to_seconds( st1.time_total ), to_seconds( st2.time_total ), cec1, cec2 );
   // break;
 }

 exp.save();
 exp.table();

 return 0;
}