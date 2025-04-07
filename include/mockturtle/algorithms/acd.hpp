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
 \file acd.hpp
 \brief Ashenhurst-Curtis decomposition

 \author Alessandro Tempia Calvino
*/

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "kitty/detail/constants.hpp"
#include "kitty/constructors.hpp"
#include "kitty/dynamic_truth_table.hpp"
#include "kitty/operations.hpp"
#include "kitty/operators.hpp"
#include "kitty/print.hpp"
#include "kitty/static_truth_table.hpp"

namespace mockturtle
{

/*! \brief Parameters for acd */
struct acd_params
{
 /*! \brief LUT size for decomposition (3 < num < 7). */
 uint32_t lut_size{ 6 };

 /*! \brief Maximum size of the free set (1 < num < 6). */
 uint32_t max_free_set_vars{ 5 };

 /*! \brief Maximum multiplicity value allowed (2 <= num < 16). */
 uint32_t max_multiplicity{ 16 };

 /*! \brief Perform only support reducing (2-level) decompositions. */
 bool support_reducing_only{ true };

 /*! \brief Use the first feasible decomposition found. */
 bool use_first{ false };

 /*! \brief If decomposition with delay profile fails, try without. */
 bool try_no_late_arrival{ false };
};

/*! \brief Statistics for acd */
struct acd_stats
{
 uint32_t num_luts{ 0 };
 uint32_t num_edges{ 0 };
 uint32_t num_levels{ 0 };
};

struct acd_result
{
 kitty::dynamic_truth_table tt;
 std::vector<uint32_t> support;
};

class acd_impl
{
private:
 struct encoding_column
 {
   uint64_t column[2];
   uint32_t cost;
   uint32_t index;
   float sort_cost;
 };

private:
 static constexpr uint32_t max_num_vars = 11;
 using STT = kitty::static_truth_table<max_num_vars>;
 using word = uint64_t;

public:
 explicit acd_impl( uint32_t num_vars, acd_params const& ps, acd_stats* pst = nullptr )
     : num_vars( num_vars ), ps( ps ), pst( pst )
 {
   std::iota( permutations.begin(), permutations.end(), 0 );
 }

 /*! \brief Runs ACD using late arriving variables */
 int run( word* ptt, unsigned delay_profile )
 {
   /* truth table is too large for the settings */
   if ( num_vars > max_num_vars )
   {
     return -1;
   }

   uint32_t late_arriving = __builtin_popcount( delay_profile );

   /* relax maximum number of free set variables if a function has more variables */
   if ( num_vars > ps.max_free_set_vars + ps.lut_size )
   {
     ps.max_free_set_vars = num_vars - ps.lut_size;
   }
   if ( late_arriving > ps.max_free_set_vars )
   {
     return -1;
     // ps.max_free_set_vars = late_arriving;
   }

   /* return a high cost if too many late arriving variables */
   if ( late_arriving > ps.lut_size - 1 || late_arriving > ps.max_free_set_vars )
   {
     return -1;
   }

   /* convert to static TT */
   init_truth_table( ptt );

   /* permute late arriving variables to be the least significant */
   reposition_late_arriving_variables( delay_profile, late_arriving );

   /* run ACD trying different bound sets and free sets */
   if ( !find_decomposition( delay_profile, late_arriving ) )
   {
     return -1;
   }

   /* return number of levels */
   return delay_profile == 0 ? 2 : 1;
 }

 int compute_decomposition()
 {
   if ( best_multiplicity == UINT32_MAX )
     return -1;

   /* compute isets */
   std::vector<STT> isets = compute_isets();

   generate_support_minimization_encodings();

   /* solves exactly only for small multiplicities */
   if ( best_multiplicity <= 4u )
     solve_min_support_exact( isets );
   else
     solve_min_support_heuristic( isets );

   /* unfeasible decomposition */
   assert( !best_bound_sets.empty() );

   return 0;
 }

 unsigned get_profile()
 {
   unsigned profile = 0;

   if ( best_free_set > num_vars )
     return -1;

   for ( uint32_t i = 0; i < best_free_set; ++i )
   {
     profile |= 1 << permutations[i];
   }

   return profile;
 }

 void get_decomposition( unsigned char* decompArray )
 {
   if ( best_free_set > num_vars )
     return;

   generate_decomposition();
   get_decomposition_abc( decompArray );
 }

 /*std::vector<acd_result> const& get_acd_result()
 {
   if ( best_free_set > num_vars )
     return;

   generate_decomposition();
   return dec_result;
 }*/

private:
 bool find_decomposition( unsigned& delay_profile, uint32_t late_arriving )
 {
   best_multiplicity = UINT32_MAX;
   best_free_set = UINT32_MAX;
   uint32_t best_cost = UINT32_MAX;
   uint32_t offset = static_cast<uint32_t>( late_arriving );
   uint32_t start = std::max( offset, 1u );

   /* perform only support reducing decomposition */
   if ( ps.support_reducing_only )
   {
     start = std::max( start, num_vars - ps.lut_size );
   }

   /* array of functions to compute the column multiplicity */
   std::function<uint32_t( STT const& tt )> column_multiplicity_fn[5] = {
       [this]( STT const& tt ) { return column_multiplicity<1u>( tt ); },
       [this]( STT const& tt ) { return column_multiplicity<2u>( tt ); },
       [this]( STT const& tt ) { return column_multiplicity<3u>( tt ); },
       [this]( STT const& tt ) { return column_multiplicity5<4u>( tt ); },
       [this]( STT const& tt ) { return column_multiplicity5<5u>( tt ); } };

   /* find a feasible AC decomposition */
   for ( uint32_t i = start; i <= ps.lut_size - 1 && i <= ps.max_free_set_vars; ++i )
   {
     auto [tt_p, perm, multiplicity] = enumerate_iset_combinations_offset( i, offset, column_multiplicity_fn[i - 1] );

     /* additional cost if not support reducing */
     uint32_t additional_cost = ( num_vars - i > ps.lut_size ) ? 128 : 0;

     /* check for feasible solution that improves the cost */
     if ( multiplicity <= ( 1 << ( ps.lut_size - i ) ) && multiplicity + additional_cost < best_cost && ( multiplicity <= ps.max_multiplicity || multiplicity == 16 ) )
     {
       best_tt = tt_p;
       permutations = perm;
       best_multiplicity = multiplicity;
       best_cost = multiplicity + additional_cost;
       best_free_set = i;

       if ( !ps.use_first && multiplicity > 2 )
       {
         continue;
       }

       break;
     }
   }

   if ( best_multiplicity == UINT32_MAX && ( !ps.try_no_late_arrival || late_arriving == 0 ) )
     return false;

   /* try without the delay profile */
   if ( best_multiplicity == UINT32_MAX )
   {
     delay_profile = 0;
     if ( ps.support_reducing_only )
     {
       start = std::max( 1u, num_vars - ps.lut_size );
     }

     for ( uint32_t i = start; i <= ps.lut_size - 1 && i <= ps.max_free_set_vars; ++i )
     {
       auto [tt_p, perm, multiplicity] = enumerate_iset_combinations_offset( i, 0, column_multiplicity_fn[i - 1] );

       /* additional cost if not support reducing */
       uint32_t additional_cost = ( num_vars - i > ps.lut_size ) ? 128 : 0;

       /* check for feasible solution that improves the cost */
       if ( multiplicity <= ( 1 << ( ps.lut_size - i ) ) && multiplicity + additional_cost < best_cost && ( multiplicity <= ps.max_multiplicity || multiplicity == 16 ) )
       {
         best_tt = tt_p;
         permutations = perm;
         best_multiplicity = multiplicity;
         best_cost = multiplicity + additional_cost;
         best_free_set = i;

         if ( !ps.use_first )
         {
           continue;
         }

         break;
       }
     }
   }

   if ( best_multiplicity == UINT32_MAX )
     return false;

   /* estimation on number of LUTs */
   if ( pst )
   {
     pst->num_luts = best_multiplicity <= 2 ? 2 : best_multiplicity <= 4 ? 3
                                              : best_multiplicity <= 8 ? 4
                                                                         : 5;
     pst->num_edges = ( pst->num_luts - 1 ) * ( num_vars - best_free_set ) + ( pst->num_luts - 1 ) + best_free_set;
   }

   return true;
 }

 void init_truth_table( word* ptt )
 {
   uint32_t const num_blocks = ( num_vars <= 6 ) ? 1 : ( 1 << ( num_vars - 6 ) );

   for ( uint32_t i = 0; i < num_blocks; ++i )
   {
     best_tt._bits[i] = ptt[i];
   }

   local_extend_to( best_tt, num_vars );
 }

 template<uint32_t free_set_size>
 uint32_t column_multiplicity( STT tt )
 {
   uint64_t multiplicity_set[4] = { 0u, 0u, 0u, 0u };
   uint32_t multiplicity = 0;
   uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
   uint64_t constexpr masks_bits[] = { 0x0, 0x3, 0xF, 0x3F };
   uint64_t constexpr masks_idx[] = { 0x0, 0x0, 0x0, 0x3 };

   /* supports up to 64 values of free set (256 for |FS| == 3)*/
   static_assert( free_set_size <= 3 );

   /* extract iset functions */
   auto it = std::begin( tt );
   for ( auto i = 0u; i < num_blocks; ++i )
   {
     for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
     {
       multiplicity_set[( *it >> 6 ) & masks_idx[free_set_size]] |= UINT64_C( 1 ) << ( *it & masks_bits[free_set_size] );
       *it >>= ( 1u << free_set_size );
     }
     ++it;
   }

   multiplicity = __builtin_popcountl( multiplicity_set[0] );

   if constexpr ( free_set_size == 3 )
   {
     multiplicity += __builtin_popcountl( multiplicity_set[1] );
     multiplicity += __builtin_popcountl( multiplicity_set[2] );
     multiplicity += __builtin_popcountl( multiplicity_set[3] );
   }

   return multiplicity;
 }

 template<uint32_t free_set_size>
 uint32_t column_multiplicity5( STT tt )
 {
   uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
   uint64_t constexpr masks[] = { 0x0, 0x3, 0xF, 0xFF, 0xFFFF, 0xFFFFFFFF };

   static_assert( free_set_size == 5 || free_set_size == 4 );

   uint32_t size = 0;
   uint64_t prev = -1;
   std::array<uint32_t, 64> multiplicity_set;

   /* extract iset functions */
   auto it = std::begin( tt );
   for ( auto i = 0u; i < num_blocks; ++i )
   {
     for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
     {
       uint32_t fs_fn = static_cast<uint32_t>( *it & masks[free_set_size] );
       if ( fs_fn != prev )
       {
         multiplicity_set[size++] = fs_fn;
         prev = fs_fn;
       }
       *it >>= ( 1u << free_set_size );
     }
     ++it;
   }

   std::sort( multiplicity_set.begin(), multiplicity_set.begin() + size );

   /* count unique */
   uint32_t multiplicity = 1;
   for ( auto i = 1u; i < size; ++i )
   {
     multiplicity += multiplicity_set[i] != multiplicity_set[i - 1] ? 1 : 0;
   }

   return multiplicity;
 }

 uint32_t column_multiplicity2( STT const& tt, uint32_t free_set_size )
 {
   assert( free_set_size <= 5 );

   uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
   uint64_t const shift = UINT64_C( 1 ) << free_set_size;
   uint64_t const mask = ( UINT64_C( 1 ) << shift ) - 1;
   uint32_t cofactors[4];
   uint32_t size = 0;

   /* extract iset functions */
   for ( auto i = 0u; i < num_blocks; ++i )
   {
     uint64_t sub = tt._bits[i];
     for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
     {
       uint32_t fs_fn = static_cast<uint32_t>( sub & mask );
       uint32_t k;
       for ( k = 0; k < size; ++k )
       {
         if ( fs_fn == cofactors[k] )
           break;
       }
       if ( k == 2 )
         return 3;
       if ( k == size )
         cofactors[size++] = fs_fn;
       sub >>= shift;
     }
   }

   return size;
 }

 inline bool combinations_offset_next( uint32_t k, uint32_t offset, uint32_t* pComb, uint32_t* pInvPerm, STT& tt )
 {
   uint32_t i;

   for ( i = k - 1; pComb[i] == num_vars - k + i; --i )
   {
     if ( i == offset )
       return false;
   }

   /* move vars */
   uint32_t var_old = pComb[i];
   uint32_t pos_new = pInvPerm[var_old + 1];
   std::swap( pInvPerm[var_old + 1], pInvPerm[var_old] );
   std::swap( pComb[i], pComb[pos_new] );
   kitty::swap_inplace( tt, i, pos_new );

   for ( uint32_t j = i + 1; j < k; j++ )
   {
     var_old = pComb[j];
     pos_new = pInvPerm[pComb[j - 1] + 1];
     std::swap( pInvPerm[pComb[j - 1] + 1], pInvPerm[var_old] );
     std::swap( pComb[j], pComb[pos_new] );
     kitty::swap_inplace( tt, j, pos_new );
   }

   return true;
 }

 template<typename Fn>
 std::tuple<STT, std::array<uint32_t, max_num_vars>, uint32_t> enumerate_iset_combinations_offset( uint32_t free_set_size, uint32_t offset, Fn&& fn )
 {
   STT tt = best_tt;

   /* TT with best cost */
   STT local_best_tt = tt;
   uint32_t best_cost = UINT32_MAX;

   assert( free_set_size >= offset );

   /* special case */
   if ( free_set_size == offset )
   {
     best_cost = fn( tt );
     return { tt, permutations, best_cost };
   }

   /* works up to 16 input truth tables */
   assert( num_vars <= 16 );

   /* Search for column multiplicity of 2 */
   if ( free_set_size == ps.lut_size - 1 )
   {
     return enumerate_iset_combinations2( free_set_size, offset );
   }

   /* bail early combination */
   uint32_t bail_early = 2;
   if ( best_multiplicity < UINT32_MAX )
     bail_early = ( best_multiplicity >> 1 ) + ( best_multiplicity & 1 );

   /* init combinations */
   uint32_t pComb[16], pInvPerm[16], bestPerm[16];
   for ( uint32_t i = 0; i < num_vars; ++i )
   {
     pComb[i] = pInvPerm[i] = i;
   }

   /* enumerate combinations */
   do
   {
     uint32_t cost = fn( tt );
     if ( cost < best_cost )
     {
       local_best_tt = tt;
       best_cost = cost;
       for ( uint32_t i = 0; i < num_vars; ++i )
       {
         bestPerm[i] = pComb[i];
       }

       if ( best_cost <= bail_early )
         break;
     }
   } while ( combinations_offset_next( free_set_size, offset, pComb, pInvPerm, tt ) );

   std::array<uint32_t, max_num_vars> res_perm;
   for ( uint32_t i = 0; i < num_vars; ++i )
   {
     res_perm[i] = permutations[bestPerm[i]];
   }

   return std::make_tuple( local_best_tt, res_perm, best_cost );
 }

 inline std::tuple<STT, std::array<uint32_t, max_num_vars>, uint32_t> enumerate_iset_combinations2( uint32_t free_set_size, uint32_t offset )
 {
   STT tt = best_tt;

   /* TT with best cost */
   STT local_best_tt = tt;
   uint32_t best_cost = ( 1 << ( ps.lut_size - free_set_size ) ) + 1;

   assert( free_set_size >= offset );

   /* init combinations */
   uint32_t pComb[16], pInvPerm[16];
   for ( uint32_t i = 0; i < num_vars; ++i )
   {
     pComb[i] = pInvPerm[i] = i;
   }

   /* enumerate combinations */
   std::array<uint32_t, max_num_vars> res_perm;

   do
   {
     uint32_t cost = column_multiplicity2( tt, free_set_size );
     if ( cost <= 2 )
     {
       local_best_tt = tt;
       best_cost = cost;
       for ( uint32_t i = 0; i < num_vars; ++i )
       {
         res_perm[i] = permutations[pComb[i]];
       }
       return std::make_tuple( local_best_tt, res_perm, best_cost );
     }
   } while ( combinations_offset_next( free_set_size, offset, pComb, pInvPerm, tt ) );

   return std::make_tuple( local_best_tt, res_perm, UINT32_MAX );
 }

 std::vector<STT> compute_isets( bool verbose = false )
 {
   /* construct isets involved in multiplicity */
   uint32_t isets_support = num_vars - best_free_set;
   std::vector<STT> isets( best_multiplicity );

   /* construct isets */
   std::unordered_map<uint64_t, uint32_t> column_to_iset;
   STT tt = best_tt;
   uint32_t offset = 0;
   uint32_t num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
   uint64_t constexpr masks[] = { 0x0, 0x3, 0xF, 0xFF, 0xFFFF, 0xFFFFFFFF };

   auto it = std::begin( tt );
   for ( auto i = 0u; i < num_blocks; ++i )
   {
     for ( auto j = 0; j < ( 64 >> best_free_set ); ++j )
     {
       uint64_t val = *it & masks[best_free_set];

       if ( auto el = column_to_iset.find( val ); el != column_to_iset.end() )
       {
         isets[el->second]._bits[i / ( 1u << best_free_set )] |= UINT64_C( 1 ) << ( j + offset );
       }
       else
       {
         isets[column_to_iset.size()]._bits[i / ( 1u << best_free_set )] |= UINT64_C( 1 ) << ( j + offset );
         column_to_iset[val] = column_to_iset.size();
       }

       *it >>= ( 1u << best_free_set );
     }

     offset = ( offset + ( 64 >> best_free_set ) ) % 64;
     ++it;
   }

   /* extend isets to cover the whole truth table */
   for ( STT& iset : isets )
   {
     local_extend_to( iset, isets_support );
   }

   /* save free_set functions */
   std::vector<STT> free_set_tts( best_multiplicity );

   for ( auto const& pair : column_to_iset )
   {
     free_set_tts[pair.second]._bits[0] = pair.first;
     local_extend_to( free_set_tts[pair.second], best_free_set );
   }

   /* print isets  and free set*/
   if ( verbose )
   {
     std::cout << "iSets\n";
     uint32_t i = 0;
     for ( auto iset : isets )
     {
       kitty::print_hex( iset );
       std::cout << " of func ";
       kitty::print_hex( free_set_tts[i++] );
       std::cout << "\n";
     }
   }

   best_free_set_tts = std::move( free_set_tts );

   return isets;
 }

 void generate_decomposition()
 {
   dec_result.clear();

   uint32_t num_edges = 0;
   for ( uint32_t i = 0; i < best_bound_sets.size(); ++i )
   {
     acd_result dec;
     auto tt = best_bound_sets[i];
     auto care = best_care_sets[i];

     /* compute and minimize support for bound set variables */
     uint32_t k = 0;
     for ( uint32_t j = 0; j < num_vars - best_free_set; ++j )
     {
       if ( !kitty::has_var( tt, care, j ) )
       {
         /* fix truth table */
         adjust_truth_table_on_dc( tt, care, j );
         continue;
       }

       if ( k < j )
       {
         kitty::swap_inplace( tt, k, j );
         kitty::swap_inplace( care, k, j );
       }
       dec.support.push_back( permutations[best_free_set + j] );
       ++k;
     }

     dec.tt = kitty::shrink_to( tt, dec.support.size() );
     dec_result.push_back( dec );
     num_edges += dec.support.size() > 1 ? dec.support.size() : 0;
   }

   /* compute the decomposition for the top-level LUT */
   compute_top_lut_decomposition();

   if ( pst )
   {
     pst->num_luts = dec_result.size();
     pst->num_edges = num_edges + dec_result.back().support.size();
   }
 }

 void compute_top_lut_decomposition()
 {
   uint32_t top_vars = best_bound_sets.size() + best_free_set;
   assert( top_vars <= ps.lut_size );

   /* extend bound set functions with free_set_size LSB vars */
   kitty::dynamic_truth_table tt( top_vars );

   /* compute support */
   dec_result.emplace_back();
   for ( uint32_t i = 0; i < best_free_set; ++i )
   {
     dec_result.back().support.push_back( permutations[i] );
   }

   /* create functions for bound set */
   std::vector<kitty::dynamic_truth_table> bound_set_vars;
   auto res_it = dec_result.begin();
   uint32_t offset = 0;
   for ( uint32_t i = 0; i < best_bound_sets.size(); ++i )
   {
     bound_set_vars.emplace_back( top_vars );
     kitty::create_nth_var( bound_set_vars[i], best_free_set + i );

     /* add bound-set variables to the support, remove buffers (shared set) */
     if ( res_it->support.size() == 1 )
     {
       dec_result.back().support.push_back( res_it->support.front() );
       /* it is a NOT */
       if ( ( res_it->tt._bits[0] & 1 ) == 1 )
       {
         bound_set_vars[i] = ~bound_set_vars[i];
       }
       dec_result.erase( res_it );
       ++offset;
     }
     else
     {
       dec_result.back().support.push_back( num_vars + i - offset );
       ++res_it;
     }
   }

   /* create composition function */
   for ( uint32_t i = 0; i < best_free_set_tts.size(); ++i )
   {
     kitty::dynamic_truth_table free_set_tt = kitty::shrink_to( best_free_set_tts[i], top_vars );

     /* find MUX assignments */
     for ( uint32_t j = 0; j < bound_set_vars.size(); ++j )
     {
       /* AND with ONSET or OFFSET */
       if ( ( ( best_iset_onset[j] >> i ) & 1 ) )
       {
         free_set_tt &= bound_set_vars[j];
       }
       else if ( ( ( best_iset_offset[j] >> i ) & 1 ) )
       {
         free_set_tt &= ~bound_set_vars[j];
       }
     }

     tt |= free_set_tt;
   }

   /* add top-level LUT to result */
   dec_result.back().tt = tt;
 }

 inline void reposition_late_arriving_variables( unsigned delay_profile, uint32_t late_arriving )
 {
   uint32_t k = 0;
   for ( uint32_t i = 0; i < late_arriving; ++i )
   {
     while ( ( ( delay_profile >> k ) & 1 ) == 0 )
       ++k;

     if ( permutations[i] == k )
     {
       ++k;
       continue;
     }

     std::swap( permutations[i], permutations[k] );
     kitty::swap_inplace( best_tt, i, k );
     ++k;
   }
 }

 template<class Iterator>
 void print_perm( Iterator begin, uint32_t free_set )
 {
   std::cout << "[";
   for ( uint32_t i = 0; i < num_vars; ++i )
   {
     if ( i == free_set )
     {
       std::cout << ", ";
     }
     std::cout << *begin << " ";
     ++begin;
   }
   std::cout << "]\n";
 }

 void generate_support_minimization_encodings()
 {
   uint32_t count = 0;

   /* enable don't cares only if not a power of 2 */
   uint32_t num_combs = 2;
   if ( __builtin_popcount( best_multiplicity ) == 1 )
   {
     uint32_t num_combs_exact[4] = { 1, 3, 35, 6435 };
     for ( uint32_t i = 0; i < 4; ++i )
     {
       if ( ( best_multiplicity >> i ) == 2u )
       {
         num_combs = num_combs_exact[i];
       }
     }
     support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
     generate_support_minimization_encodings_rec<false, true>( 0, 0, 0, count, best_multiplicity >> 1, true );
     assert( count == num_combs );
     return;
   }

   /* constraint the size of a ON-set OFF-set class for a strict encoding */
   uint32_t min_set_size = 1;
   if ( best_multiplicity <= 4 )
     min_set_size = 2;
   else if ( best_multiplicity <= 8 )
     min_set_size = 4;
   else
     min_set_size = 8;
   min_set_size = best_multiplicity - min_set_size;

   if ( best_multiplicity > 8 )
   {
     uint32_t class_sizes[13] = { 3, 3, 15, 25, 35, 35, 255, 501, 957, 1749, 3003, 4719, 6435 };
     num_combs = class_sizes[best_multiplicity - 3];
     support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
     generate_support_minimization_encodings_rec<false, false>( 0, 0, 0, count, min_set_size, true );
   }
   else
   {
     uint32_t class_sizes[13] = { 6, 3, 90, 130, 105, 35, 9330, 23436, 48708, 78474, 91377, 70785, 32175 };
     num_combs = class_sizes[best_multiplicity - 3];
     support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
     generate_support_minimization_encodings_rec<true, false>( 0, 0, 0, count, min_set_size, true );
   }

   assert( count == num_combs );
 }

 template<bool enable_dcset, bool equal_size_partition>
 void generate_support_minimization_encodings_rec( uint32_t onset, uint32_t offset, uint32_t var, uint32_t& count, uint32_t min_set_size, bool first )
 {
   if ( var == best_multiplicity )
   {
     if ( equal_size_partition )
     {
       /* sets must be equally populated */
       if ( __builtin_popcount( onset ) != __builtin_popcount( offset ) )
       {
         return;
       }
     }
     else if ( __builtin_popcount( onset ) < min_set_size || __builtin_popcount( offset ) < min_set_size )
     {
       /* ON-set and OFF-set must be populated with at least min_set_size elements */
       return;
     }

     support_minimization_encodings[count][0] = onset;
     support_minimization_encodings[count][1] = offset;
     ++count;
     return;
   }

   /* var in DCSET */
   if ( enable_dcset )
   {
     generate_support_minimization_encodings_rec<enable_dcset, equal_size_partition>( onset, offset, var + 1, count, min_set_size, first );
   }

   /* move var in ONSET */
   onset |= 1 << var;
   generate_support_minimization_encodings_rec<enable_dcset, equal_size_partition>( onset, offset, var + 1, count, min_set_size, false );
   onset &= ~( 1 << var );

   /* remove symmetries */
   if ( first )
   {
     return;
   }

   /* move var in OFFSET */
   offset |= 1 << var;
   generate_support_minimization_encodings_rec<enable_dcset, equal_size_partition>( onset, offset, var + 1, count, min_set_size, false );
   offset &= ~( 1 << var );
 }

 void solve_min_support_exact( std::vector<STT> const& isets )
 {
   std::vector<encoding_column> matrix;
   matrix.reserve( support_minimization_encodings.size() );
   best_bound_sets.clear();

   /* create covering matrix */
   if ( !create_covering_matrix<false>( isets, matrix, false ) )
   {
     return;
   }

   /* solve the covering problem */
   std::array<uint32_t, 6> solution = covering_solve_exact( matrix );

   /* check for failed decomposition */
   if ( solution[0] == UINT32_MAX )
   {
     return;
   }

   /* compute best bound sets */
   uint32_t num_luts = 1 + solution[5];
   uint32_t num_levels = 2;
   uint32_t num_edges = best_free_set + solution[5];
   uint32_t isets_support = num_vars - best_free_set;
   best_care_sets.clear();
   best_iset_onset.clear();
   best_iset_offset.clear();
   for ( uint32_t i = 0; i < solution[5]; ++i )
   {
     STT tt;
     STT care;

     const uint32_t onset = support_minimization_encodings[matrix[solution[i]].index][0];
     const uint32_t offset = support_minimization_encodings[matrix[solution[i]].index][1];
     for ( uint32_t j = 0; j < best_multiplicity; ++j )
     {
       if ( ( ( onset >> j ) & 1 ) )
       {
         tt |= isets[j];
       }
       if ( ( ( offset >> j ) & 1 ) )
       {
         care |= isets[j];
       }
     }

     care |= tt;
     num_edges += matrix[solution[i]].cost & ( ( 1 << isets_support ) - 1 );

     best_bound_sets.push_back( tt );
     best_care_sets.push_back( care );
     best_iset_onset.push_back( onset );
     best_iset_offset.push_back( offset );
   }

   if ( pst )
   {
     pst->num_luts = num_luts;
     pst->num_levels = num_levels;
     pst->num_edges = num_edges;
   }
 }

 void solve_min_support_heuristic( std::vector<STT> const& isets )
 {
   std::vector<encoding_column> matrix;
   matrix.reserve( support_minimization_encodings.size() );
   best_bound_sets.clear();

   /* create covering matrix */
   if ( !create_covering_matrix<true>( isets, matrix, true ) )
   {
     return;
   }

   /* solve the covering problem: heuristic pass + local search */
   std::array<uint32_t, 6> solution = covering_solve_heuristic( matrix );

   /* check for failed decomposition */
   if ( solution[0] == UINT32_MAX )
   {
     return;
   }

   /* improve solution with local search */
   while ( covering_improve( matrix, solution ) )
     ;

   /* try to increase the encoding length */
   // if ( covering_increase_length( matrix, solution ) )
   //   std::cout << "[i] stat\n";

   /* compute best bound sets */
   uint32_t num_luts = 1 + solution[5];
   uint32_t num_levels = 2;
   uint32_t num_edges = best_free_set + solution[5];
   uint32_t isets_support = num_vars - best_free_set;
   best_care_sets.clear();
   best_iset_onset.clear();
   best_iset_offset.clear();
   for ( uint32_t i = 0; i < solution[5]; ++i )
   {
     STT tt;
     STT care;

     const uint32_t onset = support_minimization_encodings[matrix[solution[i]].index][0];
     const uint32_t offset = support_minimization_encodings[matrix[solution[i]].index][1];
     for ( uint32_t j = 0; j < best_multiplicity; ++j )
     {
       if ( ( ( onset >> j ) & 1 ) )
       {
         tt |= isets[j];
       }
       if ( ( ( offset >> j ) & 1 ) )
       {
         care |= isets[j];
       }
     }

     care |= tt;
     num_edges += matrix[solution[i]].cost & ( ( 1 << isets_support ) - 1 );

     best_bound_sets.push_back( tt );
     best_care_sets.push_back( care );
     best_iset_onset.push_back( onset );
     best_iset_offset.push_back( offset );
   }

   if ( pst )
   {
     pst->num_luts = num_luts;
     pst->num_levels = num_levels;
     pst->num_edges = num_edges;
   }
 }

 template<bool UseHeuristic>
 bool create_covering_matrix( std::vector<STT> const& isets, std::vector<encoding_column>& matrix, bool sort )
 {
   assert( best_multiplicity <= 16 );
   uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
   uint32_t iset_support = num_vars - best_free_set;

   /* insert dichotomies */
   for ( uint32_t i = 0; i < support_minimization_encodings.size(); ++i )
   {
     uint32_t const onset = support_minimization_encodings[i][0];
     uint32_t const offset = support_minimization_encodings[i][1];

     /* compute function and distinguishable seed dichotomies */
     uint64_t column[2] = { 0, 0 };
     STT tt;
     STT care;
     uint32_t pair_pointer = 0;
     for ( uint32_t j = 0; j < best_multiplicity; ++j )
     {
       auto onset_shift = ( onset >> j );
       auto offset_shift = ( offset >> j );
       if ( ( onset_shift & 1 ) )
       {
         tt |= isets[j];
       }

       if ( ( offset_shift & 1 ) )
       {
         care |= isets[j];
       }

       /* compute included seed dichotomies */
       for ( uint32_t k = j + 1; k < best_multiplicity; ++k )
       {
         /* if are in diffent sets */
         if ( ( ( ( onset_shift & ( offset >> k ) ) | ( ( onset >> k ) & offset_shift ) ) & 1 ) )
         {
           column[pair_pointer >> 6u] |= UINT64_C( 1 ) << ( pair_pointer & 0x3F );
         }

         ++pair_pointer;
       }
     }

     care |= tt;

     /* compute cost */
     uint32_t cost = 0;
     for ( uint32_t j = 0; j < iset_support; ++j )
     {
       cost += has_var_support( tt, care, iset_support, j ) ? 1 : 0;
       // if ( !has_var_support( tt, care, iset_support, j ) )
       // {
       //   /* adjust truth table and care set */
       //   adjust_truth_table_on_dc_support( tt, care, iset_support, j );
       //   continue;
       // }
       // ++cost;
     }

     /* discard solutions with support over LUT size */
     if ( cost > ps.lut_size )
       continue;

     /* buffers have zero cost */
     if ( cost == 1 )
       cost = 0;

     float sort_cost = 0;
     if constexpr ( UseHeuristic )
     {
       sort_cost = 1.0f / ( __builtin_popcountl( column[0] ) + __builtin_popcountl( column[1] ) );
     }
     else
     {
       sort_cost = cost + ( ( combinations - __builtin_popcountl( column[0] + __builtin_popcountl( column[1] ) ) ) << num_vars );
     }

     /* insert */
     matrix.emplace_back( encoding_column{ { column[0], column[1] }, cost, i, sort_cost } );
   }

   if ( !sort )
   {
     return true;
   }

   if constexpr ( UseHeuristic )
   {
     std::sort( matrix.begin(), matrix.end(), [&]( auto const& a, auto const& b ) {
       return a.cost < b.cost;
     } );
   }
   else
   {
     std::sort( matrix.begin(), matrix.end(), [&]( auto const& a, auto const& b ) {
       return a.sort_cost < b.sort_cost;
     } );
   }

   return true;
 }

 std::array<uint32_t, 6> covering_solve_exact( std::vector<encoding_column>& matrix )
 {
   /* last value of res contains the size of the bound set */
   std::array<uint32_t, 6> res = { UINT32_MAX };
   uint32_t best_cost = UINT32_MAX;
   uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;

   assert( best_multiplicity <= 4 );

   /* determine the number of needed loops*/
   if ( best_multiplicity <= 2 )
   {
     res[5] = 1;
     res[0] = 0;
   }
   else if ( best_multiplicity <= 4 )
   {
     res[5] = 2;
     for ( uint32_t i = 0; i < matrix.size() - 1; ++i )
     {
       for ( uint32_t j = 1; j < matrix.size(); ++j )
       {
         /* filter by cost */
         if ( matrix[i].cost + matrix[j].cost >= best_cost )
           continue;

         /* check validity */
         if ( __builtin_popcountl( matrix[i].column[0] | matrix[j].column[0] ) + __builtin_popcountl( matrix[i].column[1] | matrix[j].column[1] ) == combinations )
         {
           res[0] = i;
           res[1] = j;
           best_cost = matrix[i].cost + matrix[j].cost;
         }
       }
     }
   }

   return res;
 }

 std::array<uint32_t, 6> covering_solve_heuristic( std::vector<encoding_column>& matrix )
 {
   /* last value of res contains the size of the bound set */
   std::array<uint32_t, 6> res = { UINT32_MAX };
   uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
   uint64_t column0 = 0, column1 = 0;

   uint32_t best = 0;
   float best_cost = std::numeric_limits<float>::max();
   for ( uint32_t i = 0; i < matrix.size(); ++i )
   {
     if ( matrix[i].sort_cost < best_cost )
     {
       best = i;
       best_cost = matrix[i].sort_cost;
     }
   }

   /* select */
   column0 = matrix[best].column[0];
   column1 = matrix[best].column[1];
   std::swap( matrix[0], matrix[best] );

   /* get max number of BS's */
   uint32_t iter = 1;

   while ( iter < ps.lut_size - best_free_set && __builtin_popcountl( column0 ) + __builtin_popcountl( column1 ) != combinations )
   {
     /* select column that minimizes the cost */
     best = 0;
     best_cost = std::numeric_limits<float>::max();
     for ( uint32_t i = iter; i < matrix.size(); ++i )
     {
       float local_cost = 1.0f / ( __builtin_popcountl( matrix[i].column[0] & ~column0 ) + __builtin_popcountl( matrix[i].column[1] & ~column1 ) );
       if ( local_cost < best_cost )
       {
         best = i;
         best_cost = local_cost;
       }
     }

     column0 |= matrix[best].column[0];
     column1 |= matrix[best].column[1];
     std::swap( matrix[iter], matrix[best] );
     ++iter;
   }

   if ( __builtin_popcountl( column0 ) + __builtin_popcountl( column1 ) == combinations )
   {
     for ( uint32_t i = 0; i < iter; ++i )
     {
       res[i] = i;
     }
     res[5] = iter;
   }

   return res;
 }

 bool covering_improve( std::vector<encoding_column> const& matrix, std::array<uint32_t, 6>& solution )
 {
   /* performs one iteration of local search */
   uint32_t best_cost = 0, local_cost = 0;
   uint32_t num_elements = solution[5];
   uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
   bool improved = false;

   /* compute current cost */
   for ( uint32_t i = 0; i < num_elements; ++i )
   {
     best_cost += matrix[solution[i]].cost;
   }

   uint64_t column0, column1;
   for ( uint32_t i = 0; i < num_elements; ++i )
   {
     /* remove element i */
     local_cost = 0;
     column0 = 0;
     column1 = 0;
     for ( uint32_t j = 0; j < num_elements; ++j )
     {
       if ( j == i )
         continue;
       local_cost += matrix[solution[j]].cost;
       column0 |= matrix[solution[j]].column[0];
       column1 |= matrix[solution[j]].column[1];
     }

     /* search for a better replecemnts */
     for ( uint32_t j = 0; j < matrix.size(); ++j )
     {
       if ( __builtin_popcount( column0 | matrix[j].column[0] ) + __builtin_popcount( column1 | matrix[j].column[1] ) != combinations )
         continue;
       if ( local_cost + matrix[j].cost < best_cost )
       {
         solution[i] = j;
         best_cost = local_cost + matrix[j].cost;
         improved = true;
       }
     }
   }

   return improved;
 }

 bool covering_increase_length( std::vector<encoding_column> const& matrix, std::array<uint32_t, 6>& solution )
 {
   /* count the number of solutions with cost > 0 */
   uint32_t best_cost = 0, local_cost = 0;
   uint32_t replaceable_columns = 0;
   uint32_t num_elements = solution[5];

   for ( uint32_t i = 0; i < num_elements; ++i )
   {
     replaceable_columns += !!matrix[solution[i]].cost;
   }

   /* not possible to increase the encoding length */
   if ( replaceable_columns == 1 || ps.lut_size - best_free_set - num_elements == 0 )
     return false;

   /* extract shared set BS candidates */
   std::vector<uint32_t> ss_candidates;
   for ( uint32_t i = 0; i < matrix.size(); ++i )
   {
     if ( matrix[i].cost == 0 )
       ss_candidates.push_back( i );
   }

   // if ( ss_candidates.size() < 2 )
   //   return false;

   /* find a column and try to replace it using two shared set variables */
   uint64_t column0, column1;
   uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
   bool improved = false;

   uint32_t i = 0;
   while ( replaceable_columns > 1 && ps.lut_size - best_free_set - num_elements > 0 )
   {
     /* find a column */
     while ( i < num_elements && matrix[solution[i]].cost == 0 )
       ++i;

     if ( i == num_elements )
       break;

     /* remove the column */
     column0 = 0;
     column1 = 0;
     for ( uint32_t j = 0; j < num_elements; ++j )
     {
       if ( j == i )
         continue;

       column0 |= matrix[solution[j]].column[0];
       column1 |= matrix[solution[j]].column[1];
     }

     /* find two buffers */
     if ( covering_increase_length_find_ss( matrix, solution, ss_candidates, column0, column1, i ) )
     {
       ++num_elements;
       --replaceable_columns;
       improved = true;
     }

     ++i;
   }

   return improved;
 }

 inline bool covering_increase_length_find_ss( std::vector<encoding_column> const& matrix, std::array<uint32_t, 6>& solution, std::vector<uint32_t> const& ss_candidates, uint64_t column0, uint64_t column1, uint32_t i )
 {
   uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;

   for ( uint32_t j = 0; j < ss_candidates.size() - 1; ++j )
   {
     for ( uint32_t k = j + 1; k < ss_candidates.size(); ++k )
     {
       uint64_t test_column0 = matrix[ss_candidates[j]].column[0] | matrix[ss_candidates[k]].column[0];
       uint64_t test_column1 = matrix[ss_candidates[j]].column[1] | matrix[ss_candidates[k]].column[1];

       if ( __builtin_popcount( column0 | test_column0 ) + __builtin_popcount( column1 | test_column1 ) != combinations )
         continue;

       /* perform the substitution */
       solution[i] = ss_candidates[j];
       solution[solution[5]] = ss_candidates[k];
       ++solution[5];
       return true;
     }
   }

   return false;
 }

 void adjust_truth_table_on_dc( STT& tt, STT& care, uint32_t var_index )
 {
   assert( var_index < tt.num_vars() );
   assert( tt.num_vars() == care.num_vars() );

   if ( tt.num_vars() <= 6 || var_index < 6 )
   {
     auto it_tt = std::begin( tt._bits );
     auto it_care = std::begin( care._bits );
     while ( it_tt != std::end( tt._bits ) )
     {
       uint64_t new_bits = *it_tt & *it_care;
       *it_tt = ( ( new_bits | ( new_bits >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index] ) |
                ( ( new_bits | ( new_bits << ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections[var_index] );
       *it_care = ( *it_care | ( *it_care >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index];
       *it_care = *it_care | ( *it_care << ( uint64_t( 1 ) << var_index ) );

       ++it_tt;
       ++it_care;
     }
     return;
   }

   const auto step = 1 << ( var_index - 6 );
   for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); i += 2 * step )
   {
     for ( auto j = 0; j < step; ++j )
     {
       tt._bits[i + j] = ( tt._bits[i + j] & care._bits[i + j] ) | ( tt._bits[i + j + step] & care._bits[i + j + step] );
       tt._bits[i + j + step] = tt._bits[i + j];
       care._bits[i + j] = care._bits[i + j] | care._bits[i + j + step];
       care._bits[i + j + step] = care._bits[i + j];
     }
   }
 }

 void local_extend_to( STT& tt, uint32_t real_num_vars )
 {
   if ( real_num_vars < 6 )
   {
     auto mask = *tt.begin();

     for ( auto i = real_num_vars; i < num_vars; ++i )
     {
       mask |= ( mask << ( 1 << i ) );
     }

     std::fill( tt.begin(), tt.end(), mask );
   }
   else
   {
     uint32_t num_blocks = ( 1u << ( real_num_vars - 6 ) );
     auto it = tt.begin();
     while ( it != tt.end() )
     {
       it = std::copy( tt.cbegin(), tt.cbegin() + num_blocks, it );
     }
   }
 }

 bool has_var_support( const STT& tt, const STT& care, uint32_t real_num_vars, uint8_t var_index )
 {
   assert( var_index < real_num_vars );
   assert( real_num_vars <= tt.num_vars() );
   assert( tt.num_vars() == care.num_vars() );

   const uint32_t num_blocks = real_num_vars <= 6 ? 1 : ( 1 << ( real_num_vars - 6 ) );
   if ( real_num_vars <= 6 || var_index < 6 )
   {
     auto it_tt = std::begin( tt._bits );
     auto it_care = std::begin( care._bits );
     while ( it_tt != std::begin( tt._bits ) + num_blocks )
     {
       if ( ( ( ( *it_tt >> ( uint64_t( 1 ) << var_index ) ) ^ *it_tt ) & kitty::detail::projections_neg[var_index] & ( *it_care >> ( uint64_t( 1 ) << var_index ) ) & *it_care ) != 0 )
       {
         return true;
       }
       ++it_tt;
       ++it_care;
     }

     return false;
   }

   const auto step = 1 << ( var_index - 6 );
   for ( auto i = 0u; i < num_blocks; i += 2 * step )
   {
     for ( auto j = 0; j < step; ++j )
     {
       if ( ( ( tt._bits[i + j] ^ tt._bits[i + j + step] ) & care._bits[i + j] & care._bits[i + j + step] ) != 0 )
       {
         return true;
       }
     }
   }

   return false;
 }

 void adjust_truth_table_on_dc_support( STT& tt, STT& care, uint32_t real_num_vars, uint32_t var_index )
 {
   assert( var_index < real_num_vars );
   assert( tt.num_vars() == care.num_vars() );

   const uint32_t num_blocks = real_num_vars <= 6 ? 1 : ( 1 << ( real_num_vars - 6 ) );
   if ( real_num_vars <= 6 || var_index < 6 )
   {
     auto it_tt = std::begin( tt._bits );
     auto it_care = std::begin( care._bits );
     while ( it_tt != std::begin( tt._bits ) + num_blocks )
     {
       uint64_t new_bits = *it_tt & *it_care;
       *it_tt = ( ( new_bits | ( new_bits >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index] ) |
                ( ( new_bits | ( new_bits << ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections[var_index] );
       *it_care = ( *it_care | ( *it_care >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index];
       *it_care = *it_care | ( *it_care << ( uint64_t( 1 ) << var_index ) );

       ++it_tt;
       ++it_care;
     }
     return;
   }

   const auto step = 1 << ( var_index - 6 );
   for ( auto i = 0u; i < static_cast<uint32_t>( num_blocks ); i += 2 * step )
   {
     for ( auto j = 0; j < step; ++j )
     {
       tt._bits[i + j] = ( tt._bits[i + j] & care._bits[i + j] ) | ( tt._bits[i + j + step] & care._bits[i + j + step] );
       tt._bits[i + j + step] = tt._bits[i + j];
       care._bits[i + j] = care._bits[i + j] | care._bits[i + j + step];
       care._bits[i + j + step] = care._bits[i + j];
     }
   }
 }

 /* Decomposition format for ABC
  *
  * The record is an array of unsigned chars where:
  *   - the first unsigned char entry stores the number of unsigned chars in the record
  *   - the second entry stores the number of LUTs
  * After this, several sub-records follow, each representing one LUT as follows:
  *   - an unsigned char entry listing the number of fanins
  *   - a list of fanins, from the LSB to the MSB of the truth table. The N inputs of the original function
  *     have indexes from 0 to N-1, followed by the internal signals in a topological order
  *   - the LUT truth table occupying 2^(M-3) bytes, where M is the fanin count of the LUT, from the LSB to the MSB.
  *     A 2-input LUT, which takes 4 bits, should be stretched to occupy 8 bits (one unsigned char)
  *     A 0- or 1-input LUT can be represented similarly but it is not expected that such LUTs will be represented
  */
 void get_decomposition_abc( unsigned char* decompArray )
 {
   unsigned char* pArray = decompArray;
   unsigned char bytes = 2;

   /* write number of LUTs */
   pArray++;
   *pArray = dec_result.size();
   pArray++;

   /* write LUTs */
   for ( acd_result const& lut : dec_result )
   {
     /* write fanin size*/
     *pArray = lut.support.size();
     pArray++;
     ++bytes;

     /* write support */
     for ( uint32_t i : lut.support )
     {
       *pArray = (unsigned char)i;
       pArray++;
       ++bytes;
     }

     /* write truth table */
     uint32_t tt_num_bytes = ( lut.tt.num_vars() <= 3 ) ? 1 : ( 1 << ( lut.tt.num_vars() - 3 ) );
     tt_num_bytes = std::min( tt_num_bytes, 8u );
     for ( uint32_t i = 0; i < lut.tt.num_blocks(); ++i )
     {
       for ( uint32_t j = 0; j < tt_num_bytes; ++j )
       {
         *pArray = (unsigned char)( ( lut.tt._bits[i] >> ( 8 * j ) ) & 0xFF );
         pArray++;
         ++bytes;
       }
     }
   }

   /* write numBytes */
   *decompArray = bytes;
 }

private:
 uint32_t best_multiplicity{ UINT32_MAX };
 uint32_t best_free_set{ UINT32_MAX };
 STT best_tt;
 std::vector<STT> best_bound_sets;
 std::vector<STT> best_care_sets;
 std::vector<STT> best_free_set_tts;
 std::vector<uint64_t> best_iset_onset;
 std::vector<uint64_t> best_iset_offset;
 std::vector<acd_result> dec_result;

 std::vector<std::array<uint32_t, 2>> support_minimization_encodings;

 uint32_t num_vars;
 acd_params ps;
 acd_stats* pst;
 std::array<uint32_t, max_num_vars> permutations;
};

} // namespace mockturtle