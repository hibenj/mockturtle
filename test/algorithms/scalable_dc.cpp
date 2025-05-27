#include <catch.hpp>

#include <filesystem>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/algorithms/reconv_cut.hpp>
#include <mockturtle/algorithms/scalable_dc.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/window_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/topo_view.hpp>
#include <mockturtle/views/window_view.hpp>

using namespace mockturtle;

std::tuple<mockturtle::aig_network, std::vector<mockturtle::aig_network::node>, uint32_t> create_cut()
{
  aig_network aig{};
  const auto pi1 = aig.create_pi();
  const auto pi2 = aig.create_pi();
  const auto pi3 = aig.create_pi();
  const auto pi4 = aig.create_pi();
  const auto pi5 = aig.create_pi();
  const auto pi6 = aig.create_pi();
  const auto pi7 = aig.create_pi();
  const auto pi8 = aig.create_pi();
  const auto pi9 = aig.create_pi();
  const auto pi10 = aig.create_pi();
  const auto pi11 = aig.create_pi();
  const auto pi12 = aig.create_pi();
  const auto pi13 = aig.create_pi();
  const auto pi14 = aig.create_pi();
  //
  const auto n16 = aig.create_and( pi1, pi2 );
  const auto n15 = aig.create_and( !pi1, !pi2 );
  const auto n27 = aig.create_and( pi13, pi14 );
  const auto n25 = aig.create_and( pi11, pi12 );
  const auto n24 = aig.create_and( !pi11, !pi12 );
  const auto n30 = aig.create_and( pi9, pi10 );
  const auto n23 = aig.create_and( !pi9, !pi10 );
  const auto n34 = aig.create_and( pi7, pi8 );
  const auto n22 = aig.create_and( !pi7, !pi8 );
  const auto n20 = aig.create_and( pi5, pi6 );
  const auto n19 = aig.create_and( !pi5, !pi6 );
  const auto n40 = aig.create_and( pi3, pi4 );
  const auto n18 = aig.create_and( !pi3, !pi4 );
  //
  const auto n17 = aig.create_and( !n16, !n15 );
  const auto n26 = aig.create_and( !n25, !n24 );
  const auto n31 = aig.create_and( !n30, !n23 );
  const auto n35 = aig.create_and( !n34, !n22 );
  const auto n21 = aig.create_and( !n20, !n19 );
  const auto n41 = aig.create_and( !n40, !n18 );
  //
  const auto n28 = aig.create_and( !n27, n26 );
  const auto n29 = aig.create_and( !n28, !n24 );
  const auto n32 = aig.create_and( !n29, n31 );
  const auto n33 = aig.create_and( !n32, !n23 );
  const auto n36 = aig.create_and( !n33, n35 );
  const auto n37 = aig.create_and( !n36, !n22 );
  const auto n38 = aig.create_and( !n37, n21 );
  const auto n39 = aig.create_and( !n38, !n19 );
  const auto n42 = aig.create_and( !n39, n41 );
  const auto n43 = aig.create_and( !n42, !n18 );
  //
  const auto n44 = aig.create_and( n17, !n43 );
  const auto n45 = aig.create_and( !n44, !n15 );

  aig.create_po( n45 );

  std::vector<mockturtle::node<aig_network>> leaves{};
  leaves.push_back( aig.get_node( n15 ) );
  leaves.push_back( aig.get_node( n17 ) );
  leaves.push_back( aig.get_node( n18 ) );
  leaves.push_back( aig.get_node( n40 ) );
  leaves.push_back( aig.get_node( pi5 ) );
  leaves.push_back( aig.get_node( pi6 ) );
  leaves.push_back( aig.get_node( n37 ) );

  const auto index = aig.node_to_index( aig.get_node( n45 ) );

  return std::make_tuple( aig, leaves, index );
}

std::tuple<mockturtle::aig_network, std::vector<mockturtle::aig_network::node>, uint32_t> create_cut2()
{
  aig_network aig{};
  const auto pi1 = aig.create_pi();
  const auto pi2 = aig.create_pi();
  const auto pi3 = aig.create_pi();
  const auto pi4 = aig.create_pi();
  const auto pi5 = aig.create_pi();
  const auto n1 = aig.create_and( pi1, pi2 );
  const auto n2 = aig.create_and( !pi2, !pi3 );
  const auto n3 = aig.create_and( n1, n2 );
  const auto n4 = aig.create_and( pi4, pi5 );
  const auto n5 = aig.create_and( n3, n4 );
  aig.create_po( n5 );

  std::vector<mockturtle::node<aig_network>> leaves{};
  leaves.push_back( aig.get_node( n1 ) );
  leaves.push_back( aig.get_node( n2 ) );
  leaves.push_back( aig.get_node( n4 ) );

  const auto index = aig.node_to_index( aig.get_node( n5 ) );

  return std::make_tuple( aig, leaves, index );
}

std::tuple<mockturtle::aig_network, std::vector<mockturtle::aig_network::node>, uint32_t> create_cut3()
{
  aig_network aig{};
  const auto pi1 = aig.create_pi();
  const auto pi2 = aig.create_pi();
  const auto pi3 = aig.create_pi();
  const auto pi4 = aig.create_pi();
  const auto pi5 = aig.create_pi();
  const auto pi6 = aig.create_pi();
  const auto pi7 = aig.create_pi();
  const auto pi8 = aig.create_pi();
  const auto pi9 = aig.create_pi();
  const auto pi10 = aig.create_pi();
  const auto pi11 = aig.create_pi();
  const auto pi12 = aig.create_pi();
  const auto pi13 = aig.create_pi();
  const auto pi14 = aig.create_pi();
  //
  const auto n29 = aig.create_and(!pi11, !pi12);
  const auto n31 = aig.create_and(!pi14, !pi8);
  const auto n17 = aig.create_and(!pi5, !pi8);
  const auto n25 = aig.create_and(pi5, !pi4);
  const auto n18 = aig.create_and(pi4, !pi9);
  const auto n15 = aig.create_and(pi3, !pi7);
  const auto n21 = aig.create_and(pi9, !pi6);
  //
  const auto n30 = aig.create_and(!n29, !pi4);
  const auto n32 = aig.create_and(pi13, !n31);
  const auto n19 = aig.create_and(n17, n18);
  const auto n16 = aig.create_and(n15, pi6);
  //
  const auto n33 = aig.create_and(!n32, !pi6);
  const auto n20 = aig.create_and(n19, n16);
  const auto n22 = aig.create_and(!n16, !n21);
  //
  const auto n34 = aig.create_and(!n30, !n33);
  const auto n23 = aig.create_and(pi8, !n22);
  //
  const auto n35 = aig.create_and(!n34, !pi3);
  const auto n24 = aig.create_and(!n23, !pi10);
  //
  const auto n26 = aig.create_and(n25, !n24);
  const auto n27 = aig.create_and(!n26, !n20);
  const auto n28 = aig.create_and(!pi2, !n27);
  const auto n36 = aig.create_and(!n35, !n28);
  const auto n37 = aig.create_and(!pi1, !n36);

  aig.create_po(n37);

  std::vector<mockturtle::node<aig_network>> leaves{};
  leaves.push_back( aig.get_node( pi1 ) );
  leaves.push_back( aig.get_node( pi2 ) );
  leaves.push_back( aig.get_node( n27 ) );
  leaves.push_back( aig.get_node( pi3 ) );
  leaves.push_back( aig.get_node( n33 ) );
  leaves.push_back( aig.get_node( n29 ) );
  leaves.push_back( aig.get_node( pi4 ) );

  const auto index = aig.node_to_index( aig.get_node( n37 ) );

  /*std::ostringstream out;
  write_dot( aig, out );
  std::cout << out.str() << std::endl;*/

  return std::make_tuple( aig, leaves, index );
}

std::tuple<mockturtle::aig_network, std::vector<mockturtle::aig_network::node>, uint32_t> create_cut_free()
{
  aig_network aig{};
  const auto pi1 = aig.create_pi();
  const auto pi2 = aig.create_pi();
  const auto n1 = aig.create_and( pi1, pi2 );
  const auto n2 = aig.create_and( !pi1, !pi2 );
  const auto n3 = aig.create_and( n1, n2 );
  aig.create_po( n3 );

  std::vector<mockturtle::node<aig_network>> leaves{};
  leaves.push_back( aig.get_node( n1 ) );
  leaves.push_back( aig.get_node( n2 ) );

  const auto index = aig.node_to_index( aig.get_node( n3 ) );

  return std::make_tuple( aig, leaves, index );
}

uint32_t count_zeros( const kitty::dynamic_truth_table& tt, uint32_t cut_size )
{
  uint32_t const num_blocks = ( cut_size > 6 ) ? ( 1u << ( cut_size - 6 ) ) : 1;
  uint32_t count = 0;
  for ( uint32_t i = 0; i < num_blocks; ++i )
  {
    const auto c = tt._bits[i];
    count += static_cast<uint32_t>( 64 - std::__popcount( c ) ); // count unset bits
  }

  return count;
}

template<typename Ntk>
kitty::dynamic_truth_table simulate_window( Ntk ntk, uint32_t index, const std::vector<mockturtle::node<Ntk>>& leaves )
{
  static constexpr uint32_t window_size = 14;
  static constexpr uint32_t max_window_size = 14;

  reconvergence_driven_cut_parameters rps;
  rps.max_leaves = window_size; // ps.window_size;
  reconvergence_driven_cut_statistics rst;
  detail::reconvergence_driven_cut_impl<Ntk, false, false> reconv_cuts( ntk, rps, rst );

  mockturtle::color_view<Ntk> color_ntk{ ntk };

  std::vector<mockturtle::node<Ntk>> roots = { index };
  auto const extended_leaves = reconv_cuts.run( roots ).first;

 /* std::vector<mockturtle::node<Ntk>> test_leaves{};
  test_leaves.push_back( 1 );
  test_leaves.push_back( 2 );
  test_leaves.push_back( 3 );
  test_leaves.push_back( 4 );
  test_leaves.push_back( 5 );
  test_leaves.push_back( 6 );
  test_leaves.push_back( 39 );*/

  /*auto depth_ntk = mockturtle::depth_view( mockturtle::fanout_view( color_ntk ) );
  create_window_impl windowing( depth_ntk );
  const auto res = windowing.run( roots[0], window_size, 8 );*/

  std::vector<mockturtle::node<Ntk>> gates{ collect_nodes( color_ntk, extended_leaves, roots ) };
  window_view window_ntk{ color_ntk, extended_leaves, roots, gates };

  /*std::ostringstream out;
  write_dot( window_ntk, out );
  std::cout << out.str() << std::endl;

  std::vector<mockturtle::node<Ntk>> gates2{ collect_nodes( color_ntk, leaves, roots ) };
  window_view window_ntk2{ color_ntk, leaves, roots, gates2 };

  std::cout << "\n";
  std::cout << "\n";
  std::cout << "\n";
  std::ostringstream out2;
  write_dot( window_ntk2, out2 );
  std::cout << out2.str() << std::endl;*/

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

template<typename Ntk>
void print_fi_relations( Ntk const& ntk )
{
  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_constant( n ) )
      return; // skip constant node

    std::cout << "Node " << ntk.node_to_index( n ) << " <- { ";

    ntk.foreach_fanin( n, [&]( auto const& fi ) {
      auto fn = ntk.get_node( fi );
      std::cout << ntk.node_to_index( fn ) << " ";
    } );

    std::cout << "}\n";
  } );
}

TEST_CASE( "Find single window", "[scalable_dc]" )
{
  auto start = std::chrono::high_resolution_clock::now();
  auto [aig, leaves, index] = create_cut3();

  static constexpr unsigned cut_size = 7u;
  cut_enumeration_params ps_c;
  ps_c.cut_size = cut_size;
  ps_c.cut_limit = 40u;
  cut_enumeration_stats st_c;
  const auto aig_topo = mockturtle::topo_view( aig );
  const auto cuts = fast_cut_enumeration<decltype( aig_topo ), cut_size, true, cut_enumeration_params>( aig_topo, ps_c, &st_c );
  for ( auto& cut : cuts.cuts( index ) )
  {
    if ( cut->size() != cut_size )
    {
      // Ignore cuts which are not maximum size
      continue;
    }
    std::vector<mockturtle::node<aig_network>> _leaves;
    for ( const auto& l : *cut )
    {
      _leaves.push_back( l );
    }
    const auto care = simulate_window( aig, index, _leaves );
    const auto scalable_care = mockturtle::scalable_dc(aig, _leaves);

    const auto zeros = count_zeros( care, _leaves.size() );
    //std::cout << "Num Zeros: " << zeros << std::endl;
    const auto zeros_s = count_zeros( scalable_care, _leaves.size() );
    //std::cout << "Num scalable Zeros: " << zeros_s << std::endl;

    if ( zeros > zeros_s )
    {
      std::cout << "Error\n";
    }
  }

  const auto care = simulate_window( aig, index, leaves );

  const auto zeros = count_zeros( care, leaves.size() );
  std::cout << "Num Zeros: " << zeros << std::endl;
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Execution time: " << elapsed_seconds.count() << "s\n";
}

TEST_CASE( "Find multiple windows", "[scalable_dc]" )
{
  auto start = std::chrono::high_resolution_clock::now();
  auto [aig, leaves, index] = create_cut3();

  print_fi_relations( aig );

  const auto care = mockturtle::scalable_dc(aig, leaves);

  const auto zeros = count_zeros( care, leaves.size() );
  std::cout << "Num Zeros: " << zeros << std::endl;
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Execution time: " << elapsed_seconds.count() << "s\n";
}