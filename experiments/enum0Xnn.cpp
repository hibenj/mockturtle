//
// Created by benjamin on 4/17/24.
//

#ifndef MOCKTURTLE_PLANAR_FOUR_INPUT_HPP
#define MOCKTURTLE_PLANAR_FOUR_INPUT_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/npn.hpp>
#include <unordered_set>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>

// define library for input gates CHECK

// read in library to catch all available gates CHECK
// the AND and INV have to be included
// for all 2-input Boolean functions there is a planar solution iff the XOR gate is also included

// the algorithm could be generalized by adding all 2-input Boolean functions to the library

// define combination function
/*rules and restrictions of combinations:
 * 1. for n inputs the next level has max ( n - 1 ) 2-input functions ( AND, OR, XOR)
 * 2. Inverters should be implemented as inverting edges
 * 3. the ordering of inputs impacts the ordering of nodes in the next level -> restrictions to the fan-in/fan-ut relations
 * */

// how many layers are reasonable until no new Boolean function can be found?

#include <bitset>
#include <kitty/kitty.hpp>
#include <variant>
#include <vector>

using TT = kitty::dynamic_truth_table;
using tt_hash = kitty::hash<TT>;

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

class EnumNtk
{
public:
  using signal = uint64_t;
  using TruthTable = kitty::dynamic_truth_table;
  using binary_t = std::function<TruthTable(const TruthTable&, const TruthTable&)>;
  using ternary_t = std::function<TruthTable(const TruthTable&, const TruthTable&, const TruthTable&)>;
  using operation_t = std::variant<binary_t, ternary_t>;

  struct Pi : public TruthTable {
    TruthTable tt; // Each Pi has a unique truth table
    signal index;
    const uint64_t level = 0;
  };

  struct Edge {
    signal source;
    signal dest;
    bool inv;
  };

  struct Node {
    std::array<signal, 3> fanin_node = {0, 0, 0}; // default initialize all fanins to 0
    std::array<signal, 3> fanin_edge = {0, 0, 0};
    bool is_binary = true; // a node is binary by default
    signal index;
    uint64_t level;
    uint64_t rank;
    operation_t operation;
    TruthTable output_function; // the node's output
  };

  std::vector<Pi> pis;  // primary inputs
  std::vector<Node> nodes;  // nodes
  std::vector<Edge> edges;  // edges
  std::vector<uint64_t> nodes_per_level_count;

  // Create Primary Input
  void create_pis(int K) {
    for(int i = 0; i < K; ++i) {
      TruthTable tt(K);
      kitty::create_nth_var(tt, i);

      Pi new_pi;
      new_pi.tt = tt;
      new_pi.index = static_cast<signal>(pis.size());
      pis.push_back(new_pi);
    }
  }

  TruthTable GetTruthTableOfFanin(signal fanin) {
    if (is_pi(fanin)) {
      return pi_at(fanin).tt;
    }
    return node_at(fanin).output_function;
  }

  TruthTable get_ntk_fct() {
    if (!nodes.empty()) {
      TruthTable output_function = nodes.back().output_function;  // access last element's output_function
      return output_function;
    }
  }

  signal create_edge( signal source, signal dest, bool inv = false ) {
    Edge new_edge;
    new_edge.source = source;
    new_edge.dest = dest;
    new_edge.inv = inv;
    edges.push_back(new_edge);
    return edges.size() - 1; // return index of last element
  }

  // Method to update node properties
  void update_node_properties(Node& node, operation_t& operation) {
    node.operation = operation; // Move operation into the new node
    if (!pis.empty()) {
      node.index = static_cast<signal>(pis.size() + nodes.size()); // Index starts after the indexes of pis
    } else {
      node.index = static_cast<signal>(nodes.size()); // The current size of the nodes vector is the new node's index
    }

    // Calculate the level of the new node
    signal max_level = 0;
    for(int i = 0; i < (node.is_binary ? 2 : 3); ++i) {
      signal level = get_level(node.fanin_node[i]);
      if(level > max_level) {
        max_level = level;
      }
    }

    node.level = max_level + 1; // The new node's level is the maximum level of its fanins plus one

    // If needed, increase the number of levels tracked
    if(nodes_per_level_count.size() <= node.level) {
      nodes_per_level_count.resize(node.level + 1, 0);
    }

    // Assign the new node's rank and increment the count for this level
    node.rank = nodes_per_level_count[node.level]++;
  }

  void initialize_node(Node& new_node, operation_t& operation, const signal* fanin_nodes, const size_t fanin_count)
  {
    for(size_t i = 0; i < fanin_count; ++i) {
      new_node.fanin_node[i] = fanin_nodes[i];
    }
    new_node.is_binary = (fanin_count == 2);
    update_node_properties(new_node, operation);
    for(size_t i = 0; i < fanin_count; ++i) {
      new_node.fanin_edge[i] = create_edge( fanin_nodes[i], new_node.index );
    }
    nodes.push_back(new_node); // Add new node into nodes vector
  }

  signal create_node(operation_t& operation, const signal fi0, const signal fi1)
  {
    Node new_node;
    const signal fanin_nodes[2] = { fi0, fi1 };
    std::visit(overloaded {
                    [this, &new_node, &fanin_nodes](binary_t& func) {
                      new_node.output_function = func(GetTruthTableOfFanin(fanin_nodes[0]), GetTruthTableOfFanin(fanin_nodes[1]));
                      return;
                    },
                    [](auto&) -> void {
                      throw std::invalid_argument("Unexpected ternary operation provided for binary node creation");
                    }
                }, operation);

    initialize_node(new_node, operation, fanin_nodes, 2);

    return new_node.index;
  }

  signal create_node(operation_t& operation, const signal fi0, const signal fi1, const signal fi2)
  {
    Node new_node;
    const signal fanin_nodes[3] = { fi0, fi1, fi2 };
    std::visit(overloaded {
                    [this, &new_node, &fanin_nodes](ternary_t& func) {
                      new_node.output_function = func(GetTruthTableOfFanin(fanin_nodes[0]), GetTruthTableOfFanin(fanin_nodes[1]), GetTruthTableOfFanin(fanin_nodes[2]));
                      return;
                    },
                    [](auto&) -> void {
                      throw std::invalid_argument("Unexpected binary operation provided for ternary node creation");
                    }
                }, operation);

    initialize_node(new_node, operation, fanin_nodes, 3);

    return nodes.size() + pis.size() - 1;
  }

  // For binary operation
  void assign_node_operation_binary(Node& node, binary_t& operation)
  {
    node.operation = operation;
    node.output_function = operation(GetTruthTableOfFanin(node.fanin_node[0]), GetTruthTableOfFanin(node.fanin_node[1]));
  }

  // For ternary operation
  void assign_node_operation_ternary(Node& node, ternary_t& operation)
  {
    node.operation = operation;
    node.output_function = operation(GetTruthTableOfFanin(node.fanin_node[0]), GetTruthTableOfFanin(node.fanin_node[1]), GetTruthTableOfFanin(node.fanin_node[2]));
  }

  Pi& pi_at(signal index) {
    if (index >= pis.size())
      throw std::out_of_range("Index out of range");
    return pis[index];
  }

  Node& node_at(signal index) {
    if (index < pis.size() || index >= (pis.size() + nodes.size()))
      throw std::out_of_range("Index out of range");
    return nodes[index - pis.size()];
  }

  Edge& edge_at(signal index) {
    if (index >= (edges.size()))
      throw std::out_of_range("Index out of range");
    return edges[index];
  }

  // Function to get the level of an index
  signal get_level(signal index) {
    if(is_pi(index)) {
      return pis[index].level;
    } else if(is_node(index)) {
      return nodes[index - pis.size()].level;
    }
    throw std::out_of_range("Index out of range");
  }

  void recompute_node_function( signal node )
  {
    Node& n = node_at(node);
      std::visit(
        overloaded {
            [this, &n]( EnumNtk::binary_t& func) {
              n.output_function = func(
              edge_at(n.fanin_edge[0]).inv ? ~GetTruthTableOfFanin(n.fanin_node[0]) : GetTruthTableOfFanin(n.fanin_node[0]),
                  edge_at(n.fanin_edge[1]).inv ? ~GetTruthTableOfFanin(n.fanin_node[1]) : GetTruthTableOfFanin(n.fanin_node[1])
                  );
              // n.output_function = func(GetTruthTableOfFanin(n.fanin_node[0]), GetTruthTableOfFanin(n.fanin_node[1]));
            },
            [this, &n]( EnumNtk::ternary_t& func) {
              n.output_function = func(
                  edge_at(n.fanin_edge[0]).inv ? ~GetTruthTableOfFanin(n.fanin_node[0]) : GetTruthTableOfFanin(n.fanin_node[0]),
                  edge_at(n.fanin_edge[1]).inv ? ~GetTruthTableOfFanin(n.fanin_node[1]) : GetTruthTableOfFanin(n.fanin_node[1]),
                  edge_at(n.fanin_edge[2]).inv ? ~GetTruthTableOfFanin(n.fanin_node[2]) : GetTruthTableOfFanin(n.fanin_node[2])
              );
            }
        },
        n.operation // This should be an instance of either binary_t or ternary_t
    );
  }

  void recompute_node_functions(){
    foreach_node([this] (const auto n)
    {
      recompute_node_function( n );
    });
  }

  template<typename Fn>
  void foreach_pi(Fn&& fn) const {
    auto r = mockturtle::range<uint64_t>( 0, pis.size() );
    mockturtle::detail::foreach_element( r.begin(), r.end(), fn );
  }

  template<typename Fn>
  void foreach_node(Fn&& fn) const {
    auto start_index = pis.size();
    auto r = mockturtle::range<uint64_t>( start_index, nodes.size() + start_index );
    mockturtle::detail::foreach_element( r.begin(), r.end(), fn );
  }

  template<typename Fn>
  void foreach_edge(Fn&& fn) const {
    auto r = mockturtle::range<uint64_t>( edges.size() );
    mockturtle::detail::foreach_element( r.begin(), r.end(), fn );
  }

  // Modify Node Function
  void modify_node_operation(signal index, operation_t&& new_operation) {
    node_at(index).operation = std::move(new_operation);
  }

  [[nodiscard]]  bool is_pi(signal index) const {
    return index < pis.size();
  }

  [[nodiscard]]  bool is_node(signal index) const {
    return index >= pis.size() && index < (pis.size() + nodes.size());
  }

  int num_fis(const Node& node) {
    return node.is_binary ? 2 : 3;
  }

  int num_pis() {
    return pis.size();
  }

  int num_nodes() {
    return nodes.size();
  }
};

int iter{0};

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
      printf("-> found ");
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
void try_function( std::unordered_set<TT, kitty::hash<TT>> & classes, TT tt )
{
  update_classes( classes, tt );
}

// Recursive function to generate combinations
void generate_combinations(std::vector<std::vector<uint64_t>>& input,
                            std::vector<std::vector<std::vector<uint64_t>>>& output,
                            std::vector<std::vector<uint64_t>> temp,
                            std::vector<uint64_t> start_nodes, int index)
{
  if(!temp.empty()) {
    output.push_back(temp);
  }
  for(size_t i = index; i < input.size(); ++i) {
    auto input_iter = std::find(start_nodes.begin(), start_nodes.end(), input[i].front());
    if (temp.empty() || input_iter >= std::find(start_nodes.begin(), start_nodes.end(), temp.back().back())){
      temp.push_back(input[i]);
      generate_combinations(input, output, temp, start_nodes, static_cast<int>(i) + 1);
      temp.pop_back();
    }
  }
}

void compute_network_structures( std::vector<EnumNtk>& Networks, std::vector<uint64_t> start_nodes, EnumNtk enum_ntk, int cur_lvl, std::vector<EnumNtk::operation_t>& operations )
{
  // Code
  ++cur_lvl;
  // Debug output
  std::cout << "Called with: " << cur_lvl << std::endl;
  printf("\nPis: ");
  for(const auto &val : start_nodes) {
    std::cout << val << " ";
  }
  printf("\n");
  // Code
  // create all three input adjacent nodes
  std::vector<std::vector<uint64_t>> adjacent_nodes;
  // Create all two and three input adjacent nodes
  for (size_t i = 0; i < start_nodes.size() - 1; ++i) {
    std::vector<uint64_t> temp_2 = { start_nodes[i], start_nodes[i + 1] };
    adjacent_nodes.push_back(temp_2);
    if( i < start_nodes.size() - 2 )
    {
      std::vector<uint64_t> temp_3 = { start_nodes[i], start_nodes[i + 1], start_nodes[i + 2] };
      adjacent_nodes.push_back(temp_3);
    }
  }
  // Debug output
  printf("Before erasing\n");
  for( const auto &innerVec : adjacent_nodes ) {
    std::cout << "(";
    for( const auto &val : innerVec ) {
      std::cout << val << " ";
    }
    std::cout << ")" << std::endl;
  }
  // Code
  adjacent_nodes.erase( std::remove_if( adjacent_nodes.begin(), adjacent_nodes.end(),
                                        [&cur_lvl, &enum_ntk]( const auto &innerVec ) {
                                          for( const auto &val : innerVec )
                                          {
                                            if (enum_ntk.is_pi( val ) ) {
                                              if ( enum_ntk.pi_at( val ).level == cur_lvl ) {
                                                return false;
                                              }
                                            } else {
                                              if ( enum_ntk.node_at( val ).level == cur_lvl ) {
                                                return false;
                                              }
                                            }
                                          }
                                          return true;
                                        } ),
                        adjacent_nodes.end() );
  // Debug output
  printf("After erasing\n");
  for( const auto &innerVec : adjacent_nodes ) {
    std::cout << "(";
    for( const auto &val : innerVec ) {
      std::cout << val << " ";
    }
    std::cout << ")" << std::endl;
  }
  // Code
  std::vector<std::vector<std::vector<uint64_t>>> output;
  generate_combinations(adjacent_nodes, output, {}, start_nodes, 0);
  // Debug output
  for(const auto& combination : output) { // For each combination
    std::cout << "( ";
    for(const auto& pair : combination) { // For each pair in the combination
      std::cout << "( ";
      for(const auto& val : pair) { // For each value in the pair
        std::cout << val << " ";
      }
      std::cout << ") ";
    }
    std::cout << ")" << std::endl;
  }
  // Code
  std::vector<EnumNtk::signal> delete_nodes;
  for( const auto& combination : output ) { // For each combination
    auto current_nodes = start_nodes;
    auto current_ntk = enum_ntk;
    for( const auto& pair : combination ) { // For each pair in the combination
      if ( pair.size() == 2 ) { // For each value in the pair
        EnumNtk::signal new_node = current_ntk.create_node( operations[1], pair[0], pair[1] );
        std::cout << "Node created with fis: " << pair[0] << pair[1] << std::endl;
        // Replace pair[0] and pair[1] with new_node in start_nodes
        auto it1 = std::find(current_nodes.begin(), current_nodes.end(), pair[0]);
        auto it2 = std::find(current_nodes.begin(), current_nodes.end(), pair[1]);

        if (it1 != current_nodes.end() && it2 != current_nodes.end()) {
          // If it1 is already in delete_nodes, remove it from delete_nodes
          auto new_end = std::remove_if(delete_nodes.begin(), delete_nodes.end(),
                                         [&](const EnumNtk::signal &node) { return node == *it1; });
          delete_nodes.erase(new_end, delete_nodes.end());
          *it1 = new_node;
          // Add it2 into delete_nodes
          delete_nodes.push_back(*it2);
        }
      }
      else if ( pair.size() == 3 )
      {
        EnumNtk::signal new_node = current_ntk.create_node( operations[0], pair[0], pair[1], pair[2] );
        std::cout << "Node created with fis: " << pair[0] << pair[1] << pair[2] << std::endl;
        // Replace pair[0] and pair[1] with new_node in start_nodes
        auto it1 = std::find(current_nodes.begin(), current_nodes.end(), pair[0]);
        auto it2 = std::find(current_nodes.begin(), current_nodes.end(), pair[1]);
        auto it3 = std::find(current_nodes.begin(), current_nodes.end(), pair[2]);

        if (it1 != current_nodes.end() && it2 != current_nodes.end()) {
          // If it1 is already in delete_nodes, remove it from delete_nodes
          auto new_end = std::remove_if(delete_nodes.begin(), delete_nodes.end(),
                                         [&](const EnumNtk::signal &node) { return node == *it1; });
          delete_nodes.erase(new_end, delete_nodes.end());
          *it1 = new_node;
          // Add it2 into delete_nodes
          delete_nodes.push_back(*it2);
          delete_nodes.push_back(*it3);
        }
      }
      else
        throw std::runtime_error( "Not allowed pair size" );
    }
    for(auto& delete_node: delete_nodes) {
      auto new_end = std::remove_if(current_nodes.begin(), current_nodes.end(),
                                     [&](const EnumNtk::signal &current_node) { return current_node == delete_node; });
      current_nodes.erase(new_end, current_nodes.end());
    }
    delete_nodes.clear();
    //Debug
    printf("\nNodes for new iteration: ");
    for(const auto &val : current_nodes) {
      std::cout << val << " ";
    }
    printf("\n");
    //Code
    if ( current_nodes.size() == 1 )
    {
      Networks.push_back( current_ntk );
      std::cout << "\nNetwork added" << std::endl;
    }
    compute_network_structures( Networks, current_nodes, current_ntk, cur_lvl, operations );
    // break;
  }
}

// Then inside your function
template<class TT>
void create_and_try_functions(std::unordered_set<TT, kitty::hash<TT>> & classes, std::vector<TT> xs,
                               int numInputs, std::vector<EnumNtk::operation_t>& operations)
{
  using namespace mockturtle;
  std::vector<EnumNtk> Networks;
  EnumNtk enum_ntk;
  enum_ntk.create_pis( numInputs );

  // Code: create start_nodes vector
  int cur_lvl = -1;
  std::vector<uint64_t> start_nodes;
  enum_ntk.foreach_pi( [&]( auto const& pi ) {
    start_nodes.push_back( pi );
  });

  compute_network_structures( Networks, start_nodes, enum_ntk, cur_lvl, operations );

  printf("Networks size %ld \n", Networks.size());

  if( !Networks.empty() )
  {
    // This is the real code
    int n = 0;
    for(auto &network: Networks)
    {
      printf("Network it %i\n", n);
      int binaryNodeCount = std::count_if(network.nodes.begin(), network.nodes.end(), [](const auto& n) { return n.is_binary; });
      int combinationCount = std::pow(2, binaryNodeCount);
      for(int combination = 0; combination < combinationCount; ++combination)
      {
        printf("Combination it %i\n", combination);
        int idx = 0;
        for(auto &n: network.nodes){  //Iterating through nodes
          if (n.is_binary){
            int operationIdx = 1 + ((combination >> idx) & 1);
            network.assign_node_operation_binary(n, std::get<0>(operations[operationIdx]));
            idx++;
          }
        }
        int edge_invs = 0;
        for (int j = 0; j < (1 << network.edges.size()); ++j)
        {
          edge_invs = j;
          size_t i = 0;
          network.foreach_edge([&, i](auto const& e) mutable {
            bool inv_value = (edge_invs & (1 << i)) != 0;
            network.edge_at(e).inv = inv_value;
            ++i;
          });
          network.recompute_node_functions();
          EnumNtk::TruthTable final_tt = network.get_ntk_fct();
          kitty::print_binary(final_tt);
          printf("\n");
          update_classes(classes, final_tt);
          // nn
          update_classes(classes, kitty::unary_not(final_tt));
        }
      }
      n++;
    }
  }
  else
    printf("No Networks created");
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
    const auto res = kitty::exact_n_canonization( tt_2 );
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
    printf("\n");
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

void process_2_inp_luts_4( std::unordered_set<TT, tt_hash>& classes )
{
  /* enumerate the p classes */
  std::unordered_set<TT, tt_hash> classes_2;
  TT tt_2(2);
  do
  {
    const auto res = kitty::exact_n_canonization( tt_2 );
    classes_2.insert( std::get<0>( res ) );
    kitty::next_inplace( tt_2 );

  } while ( !kitty::is_const0( tt_2 ) );

  printf("\n2 inp LUTs\n");
  for( auto const& tt : classes_2 )
  {
    kitty::print_binary(tt);
    kitty::dynamic_truth_table extended_tt(4);

    for (size_t i = 0; i < tt.num_bits(); ++i)
    {
      if (kitty::get_bit(tt, i))
      {
        kitty::set_bit(extended_tt, i);
        kitty::set_bit(extended_tt, i + tt.num_bits());
        kitty::set_bit(extended_tt, i + 2*tt.num_bits());
        kitty::set_bit(extended_tt, i + 3*tt.num_bits());
      }
    }
    printf("\n");
    kitty::print_binary( extended_tt );
    printf("\n");
    try_function(classes, extended_tt);
    try_function(classes, kitty::unary_not(extended_tt));
    auto tt_extended1 = permute_inputs(extended_tt, {0, 2, 1, 3});
    kitty::print_binary( tt_extended1 );
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
    printf("\n");
    tt_extended1 = permute_inputs(extended_tt, {0, 3, 1, 2});
    kitty::print_binary( tt_extended1 );
    printf("\n");
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
    printf("\n");
    tt_extended1 = permute_inputs(extended_tt, {1, 2, 0, 3});
    kitty::print_binary( tt_extended1 );
    printf("\n");
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
    printf("\n");
    tt_extended1 = permute_inputs(extended_tt, {1, 3, 0, 2});
    kitty::print_binary( tt_extended1 );
    printf("\n");
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
    printf("\n");
    tt_extended1 = permute_inputs(extended_tt, {2, 3, 0, 1});
    kitty::print_binary( tt_extended1 );
    printf("\n");
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
  }
}

void process_3_inp_luts_4( std::unordered_set<TT, tt_hash>& classes )
{
  /* enumerate the p classes */
  std::unordered_set<TT, tt_hash> classes_3;
  TT tt_3(3);
  do
  {
    const auto res = kitty::exact_n_canonization( tt_3 );
    classes_3.insert( std::get<0>( res ) );
    kitty::next_inplace( tt_3 );

  } while ( !kitty::is_const0( tt_3 ) );

  printf("\n2 inp LUTs\n");
  for( auto const& tt : classes_3 )
  {
    kitty::print_binary( tt );
    kitty::dynamic_truth_table extended_tt(4);

    for ( size_t i = 0; i < tt.num_bits(); ++i )
    {
      if ( kitty::get_bit( tt, i ) )
      {
        kitty::set_bit( extended_tt, i );
        kitty::set_bit( extended_tt, i + tt.num_bits() );
      }
    }
    try_function(classes, extended_tt);
    try_function(classes, kitty::unary_not(extended_tt));
    auto tt_extended1 = permute_inputs(extended_tt, {0, 1, 3, 2});
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));

    tt_extended1 = permute_inputs(extended_tt, {0, 2, 3, 1});
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));

    tt_extended1 = permute_inputs(extended_tt, {1, 2, 3, 0});
    try_function(classes, tt_extended1);
    try_function(classes, kitty::unary_not(tt_extended1));
  }
}

int main()
{
  static constexpr uint32_t K = 3u;

  /* enumerate the n classes */
  int f_count{0};
  std::unordered_set<TT, tt_hash> classes;
  TT tt(K);
  // This part is needed. Needs some computation time for 4 inputs though
  auto start = std::chrono::high_resolution_clock::now();
  do
  {
    // kitty::exact_n_canonization();
    const auto res = kitty::exact_n_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
    ++f_count;
  } while ( !kitty::is_const0( tt ) );
  std::size_t size_before = classes.size();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_1 = stop - start;

  // trivial cases
  kitty::dynamic_truth_table tt_zeros(K);
  update_classes(classes, tt_zeros);
  std::vector<TT> xs;
  for( int i{0}; i<K; ++i )
  {
    xs.emplace_back( K );
    kitty::create_nth_var( xs[i], i );
    update_classes(classes, xs[i]);
  }

  printf("fct_count: %i\n", f_count);
  printf("n_count: %i\n", classes.size());

  // Specify the functions used
  std::vector<EnumNtk::operation_t> operations;
  EnumNtk::binary_t binary_and_func = [](const auto& a, const auto& b){ return kitty::binary_and(a, b); };
  // EnumNtk::binary_t binary_or_func = [](const auto& a, const auto& b){ return kitty::binary_or(a, b); };
  EnumNtk::binary_t binary_xor_func = [](const auto& a, const auto& b){ return kitty::binary_xor(a, b); };
  EnumNtk::ternary_t ternary_maj_func = [](const auto& a, const auto& b, const auto& c){ return kitty::ternary_majority(a, b, c); };
  operations.push_back(ternary_maj_func);
  operations.push_back(binary_and_func);
  operations.push_back(binary_xor_func);

  start = std::chrono::high_resolution_clock::now();
  // MAIN FUNCTION
  create_and_try_functions( classes, xs, K, operations );
  stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_2 = stop - start;

  //
  process_2_inp_luts_3(classes);
  //
  // process_2_inp_luts_4(classes);
  // process_3_inp_luts_4(classes);
  std::size_t size_after = classes.size();

  printf("\n");
  printf("NOT CROSSING FREE\n");
  for( auto const& tt : classes )
  {
    kitty::print_binary( tt );
    printf("\n");
  }
  printf("\nNum Classes before: %zu", size_before);
  printf("\nNum Classes after: %zu", size_after);
  printf("\nTime taken by function1: %f", duration_1.count());
  printf("\nTime taken by function2: %f", duration_2.count());

  return 0;
}

#endif // MOCKTURTLE_PLANAR_FOUR_INPUT_HPP
