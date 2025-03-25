//
// Created by benjamin on 3/25/25.
//
#include <catch.hpp>

#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/binding_view.hpp>

#include <mockturtle/algorithms/extract_care_set_sat.hpp>

using namespace mockturtle;

TEST_CASE( "And Or example", "[extract_care_set_sat]" )
{
  aig_network aig{};
  const auto pi0 = aig.create_pi();
  const auto pi1 = aig.create_pi();
  const auto and0 = aig.create_and( pi0, pi1 );
  const auto or0 = aig.create_or( pi0, pi1 );
  const auto and1 = aig.create_and( or0, and0 );
  aig.create_po( and1 );

  static constexpr unsigned cut_size = 2u;
  cut_enumeration_params ps_c;
  ps_c.cut_size = cut_size;
  ps_c.cut_limit = 20u;
  cut_enumeration_stats st_c;
  const auto aig_topo = mockturtle::topo_view( aig );
  const auto cuts = fast_cut_enumeration<decltype( aig_topo ), cut_size, true, cut_enumeration_params>( aig_topo, ps_c, &st_c );

  const auto index = aig.node_to_index( aig.get_node( and1 ) );
  for ( auto& cut : cuts.cuts( index ) )
  {
    if ( cut->size() != cut_size )
    {
      // Ignore cuts which are not maximum size
      continue;
    }
    const auto tt = cuts.truth_table( *cut );
    const auto care_set = extract_care_set_sat( aig, cut );

    // The truth table of this is ~or0 & and0
    CHECK( tt._bits == 2u );

    // THe pattern 11 is never satisfiable at the inputs of and1 with leaves or0 and and0
    CHECK( care_set._bits[0] == 7u );

    // only check the first cut
    break;
  }
}