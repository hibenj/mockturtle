#include <iostream>
#include <string>
#include <fstream>
#include <bitset>
#include <vector>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/npn.hpp>
#include <unordered_set>

#include <mockturtle/networks/xag.hpp>

int iter{0};
using TT = kitty::dynamic_truth_table;
using tt_hash = kitty::hash<TT>;

template<class TT>
void update_classes( std::unordered_set<TT, kitty::hash<TT>> & classes, TT tt )
{
  printf("%2d:", iter++);
  TT rep = std::get<0>( kitty::exact_n_canonization( tt ) );
  for ( auto it = classes.begin(); it != classes.end(); )
  {
    if ( kitty::equal(*it, rep) || kitty::equal(*it, ~rep) )
    {
      it = classes.erase(it);
      printf("-> found ", classes.size());
      kitty::print_binary(rep);
    }
    else
    {
      ++it;
    }
  }
  printf("\n");
}

template<class TT>
TT maj( TT const& tt1, TT const& tt2, TT const& tt3 )
{
  return (tt1&tt2) | (tt2&tt3) | (tt1&tt3);
}

template<class TT>
void try_function( std::unordered_set<TT, kitty::hash<TT>> & classes, TT tt )
{
  update_classes( classes, tt );
  update_classes( classes, kitty::unary_not(tt) );
}

kitty::dynamic_truth_table permute_inputs(kitty::dynamic_truth_table tt, std::vector<uint8_t> perm) {
  uint8_t num_vars = tt.num_vars();
  kitty::dynamic_truth_table new_tt(num_vars);
  for(size_t row = 0; row < (1u << num_vars); ++row) {
    // Generate new row index by permutation of input variables
    size_t new_row = 0;
    for(uint8_t i = 0; i < num_vars; ++i) {
      new_row |= ((row >> perm[i]) & 1) << i;
    }
    // Transfer output value to new truth table
    if(kitty::get_bit(tt, row)) {
      kitty::set_bit(new_tt, new_row);
    }
  }
  return new_tt;
}

void process_2_inp_luts_3( std::unordered_set<TT, tt_hash>& classes )
{
  /* enumerate the p classes */
  std::unordered_set<TT, tt_hash> classes_2;
  TT tt_2(2);
  do
  {
    const auto res = kitty::exact_npn_canonization( tt_2 );
    classes_2.insert( std::get<0>( res ) );
    kitty::next_inplace( tt_2 );

  } while ( !kitty::is_const0( tt_2 ) );

  printf("\n2 inp LUTs\n");
  for( auto const& tt : classes_2 )
  {
    kitty::print_binary( tt );
    kitty::dynamic_truth_table extended_tt(3);

    for ( size_t i = 0; i < tt.num_bits(); ++i )
    {
      if ( kitty::get_bit( tt, i ) )
      {
        kitty::set_bit( extended_tt, i );
        kitty::set_bit( extended_tt, i + tt.num_bits() );
      }
    }
    printf("\n");
    kitty::print_binary( extended_tt );
    printf("\n");
    try_function(classes, extended_tt);
    try_function(classes, kitty::unary_not(extended_tt));
    auto tt_extended1 = permute_inputs(extended_tt, {0, 2, 1});
    kitty::print_binary( tt_extended1 );
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
    printf("\n");
    tt_extended1 = permute_inputs(extended_tt, {2, 1, 0});
    kitty::print_binary( tt_extended1 );
    printf("\n");
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
  }
}

int main()
{
  static constexpr uint32_t K = 3u;

  std::vector<TT> xs;
  for( int i{0}; i<K; ++i )
  {
    xs.emplace_back( K );
    kitty::create_nth_var( xs[i], i );
  }

  /* enumerate the p classes */

  std::unordered_set<TT, tt_hash> classes;
  TT tt(K);
  do
  {
    const auto res = kitty::exact_n_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );

  } while ( !kitty::is_const0( tt ) );
  std::size_t size_before = classes.size();

  printf("%ld\n", classes.size());

  auto fc = (xs[0]^xs[0]);
  auto fa = (xs[0]&xs[1]);
  auto fx = (xs[2]^xs[1]);
  auto fz = (xs[1]^xs[2]);

  printf("TRIVIAL\n");

  try_function( classes, fc );
  try_function( classes, xs[0] );
  try_function( classes, fa );
  try_function( classes, fx );

  /*
     *     4           5           6           7           8          9
     *                                                                X
     *                                                              /  ~
     *     A           X           A           X           A       |    A
     *    / \         / \         / \         / \         ~ \      |  ~  ~
     *   A   \       X   \       X   \       A   \       A   \     | A    \
     *  / \   \     / \   \     / \   \     / \   \     / \   \    |/ \    \
     * x   x   x   x   x   x   x   x   x   x   x   x   x   x   x   x   x    x
     *
   */

  printf("XAIGs\n");

  auto f0 = (xs[0]&xs[1])&xs[2];
  auto f1 = (xs[0]^xs[1])^xs[2];
  auto f2 = (xs[0]^xs[1])&xs[2];
  auto f3 = (xs[0]&xs[1])^xs[2];
  auto f4 = xs[2] & ~(xs[0]&xs[1]);
  auto f5 = xs[0]^( ( xs[2] | (xs[1]&xs[0]) ) );

  try_function( classes, f0 );
  try_function( classes, f1 );
  try_function( classes, f2 );
  try_function( classes, f3 );
  try_function( classes, f4 );
  try_function( classes, f5 );


  printf("XMIGs\n");
  auto maj3 = xs[0]&xs[1] | xs[1]&xs[2] | xs[2]&xs[0];

  /*       10        11      12
     *
     *                 X       M
     *               / |      /|\
     *       M      |  M     | X \
     *     / | \    |/ | \   |/ \ \
     *    x  x  x   x  x  x  x  x  x
     *
   */

  try_function( classes, maj( xs[0], xs[1], xs[2] ) );
  try_function( classes, xs[0]^maj( xs[0], xs[1], xs[2] ) );
  try_function( classes, maj(xs[0], fx, xs[2]) );

  /* enumerate the p classes */
  process_2_inp_luts_3(classes);

  printf("\n");
  printf("NOT CROSSING FREE\n");
  for( auto const& tt : classes )
  {
    kitty::print_binary( tt );
    printf("\n");
  }
  std::size_t size_after = classes.size();
  printf("\nNum Classes before: %zu", size_before);
  printf("\nNum Classes after: %zu", size_after);

  return 0;
}

