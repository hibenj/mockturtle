//
// Created by benjamin on 4/17/24.
//

#ifndef MOCKTURTLE_PLANAR_FOUR_INPUT_HPP
#define MOCKTURTLE_PLANAR_FOUR_INPUT_HPP

#include <iostream>
#include <string>
#include <fstream>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>
#include <kitty/npn.hpp>
#include <unordered_set>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>

// define library for input gates

// read in library to catch all available gates
// the AND and INV have to be included
// for all 2-input Boolean functions there is a planar solution iff the XOR gate is also included

// the algorithm could be generalized by adding all 2-input Boolean functions to the library

// define combination function
/*rules and restrictions of combinations:
 * 1. for n inputs the next level has max ( n - 1 ) 2-input functions ( AND, OR, XOR)
 * 2. for n inputs the next level has max n nodes ( 2-input functions + INVs )
 * 3. the ordering of inputs impacts the ordering of ndoes in the next level -> restrictions to the fan-in/fan-ut relations
 * */

// how many layers are reasonable until no new Boolean function can be found?

#include <kitty/kitty.hpp>
#include <variant>
#include <vector>

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

class EnumerationNetwork {
public:
  using signal = uint64_t;
  using TruthTable = kitty::dynamic_truth_table;
  using binary_t = std::function<TruthTable(const TruthTable&, const TruthTable&)>;
  using ternary_t = std::function<TruthTable(const TruthTable&, const TruthTable&, const TruthTable&)>;
  using operation_t = std::variant<binary_t, ternary_t>;

  struct Node {
    std::array<signal, 3> fanin = {0, 0, 0}; // default initialize all fanins to 0
    bool is_binary = true; // a node is binary by default
    signal index;
    uint64_t level;
    uint64_t rank;
    operation_t operation;
    TruthTable output_function; // the node's output
  };

  struct Pi : public TruthTable {
    TruthTable tt; // Each Pi has a unique truth table
    signal index;
    const uint64_t level = 0;
  };

  std::vector<Pi> pis;  // primary inputs
  std::vector<Node> nodes;  // nodes
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

  // Method to update node properties
  void update_node_properties(Node& node, operation_t&& operation) {
    node.operation = std::move(operation); // Move operation into the new node
    if (!pis.empty()) {
      node.index = static_cast<signal>(pis.size() + nodes.size()); // Index starts after the indexes of pis
    } else {
      node.index = static_cast<signal>(nodes.size()); // The current size of the nodes vector is the new node's index
    }

    // Calculate the level of the new node
    signal max_level = 0;
    for(int i = 0; i < (node.is_binary ? 2 : 3); ++i) {
      signal level = get_level(node.fanin[i]);
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

    nodes.push_back(node); // Add new node into nodes vector
  }

  TruthTable GetTruthTableOfFanin(signal fanin) {
    if (is_pi(fanin)) {
      return pi_at(fanin).tt;
    }
    return node_at(fanin).output_function;
  }

  // Method to create a binary node
  void create_node(operation_t&& operation, signal fi0, signal fi1) {
    Node new_node;
    new_node.fanin[0] = fi0;
    new_node.fanin[1] = fi1;
    new_node.is_binary = true;

    std::visit(overloaded {
                    [this, &new_node, fi0, fi1](binary_t& func) {
                      new_node.output_function = func(GetTruthTableOfFanin(fi0), GetTruthTableOfFanin(fi1));
                      return; // Explicitly defined return type as void
                    },
                    [](auto&) -> void { // Changed return type as void
                      throw std::invalid_argument("Unexpected ternary operation provided for binary node creation");
                    }
                }, operation);

    update_node_properties(new_node, std::move(operation));
  }

  // Method to create a ternary node
  void create_node(operation_t&& operation, signal fi0, signal fi1, signal fi2) {
    Node new_node;
    new_node.fanin[0] = fi0;
    new_node.fanin[1] = fi1;
    new_node.fanin[2] = fi2;
    new_node.is_binary = false;

    std::visit(overloaded {
                    [this, &new_node, fi0, fi1, fi2](ternary_t& func) {
                      new_node.output_function = func(GetTruthTableOfFanin(fi0), GetTruthTableOfFanin(fi1), GetTruthTableOfFanin(fi2));
                      return;
                    },
                    [](auto&) -> void {
                      throw std::invalid_argument("Unexpected binary operation provided for ternary node creation");
                    }
                }, operation);

    update_node_properties(new_node, std::move(operation));
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

  // Function to get the level of an index
  signal get_level(signal index) {
    if(is_pi(index)) {
      return pis[index].level;
    } else if(is_node(index)) {
      return nodes[index - pis.size()].level;
    }
    throw std::out_of_range("Index out of range");
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
};

int iter{0};

template<class TT>
void update_classes( std::unordered_set<TT, kitty::hash<TT>> & classes, TT tt )
{
  printf("%2d:", iter++);
  TT rep = std::get<0>( kitty::exact_npn_canonization( tt ) );
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

// Then inside your function
template<class TT>
void create_and_try_functions(std::unordered_set<TT, kitty::hash<TT>> & classes, std::vector<TT> xs,
                               int numInputs, std::vector<EnumerationNetwork::operation_t>& operations)
{
  using namespace mockturtle;
  EnumerationNetwork enum_ntk;
  enum_ntk.create_pis( numInputs );

  std::visit(
      overloaded {
          [&enum_ntk](EnumerationNetwork::binary_t& func) {
            enum_ntk.create_node(std::move(func), 0, 1);
          },
          [&enum_ntk](EnumerationNetwork::ternary_t& func) {
            enum_ntk.create_node(std::move(func), 0, 1, 2);
          }
      },
      operations[0] // This should be an instance of either binary_t or ternary_t
  );
  int numNodes = enum_ntk.nodes.size();
  kitty::dynamic_truth_table xor_o(2);
  create_from_binary_string( xor_o, "0110" );
  kitty::print_binary( xor_o );
  printf("\n1: ");
  kitty::print_binary( enum_ntk.pi_at(0).tt );
  printf("\n2: ");
  kitty::print_binary( enum_ntk.pi_at(1).tt );
  printf("\n3: ");
  kitty::print_binary( enum_ntk.node_at(2).output_function );
}

int main()
{
  using TT = kitty::dynamic_truth_table;
  static constexpr uint32_t K = 2u;

  std::vector<TT> xs;
  for( int i{0}; i<K; ++i )
  {
    xs.emplace_back( K );
    kitty::create_nth_var( xs[i], i );
  }

  /* enumerate the npn classes */
  int f_count{0};
  using tt_hash = kitty::hash<TT>;
  std::unordered_set<TT, tt_hash> classes;
  TT tt(K);
  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
    ++f_count;
  } while ( !kitty::is_const0( tt ) );

  printf("fct_count: %i\n", f_count);
  printf("npn_count: %i\n", classes.size());

  // Specify the functions used
  std::vector<EnumerationNetwork::operation_t> operations;
  EnumerationNetwork::binary_t binary_and_func = [](const auto& a, const auto& b){ return kitty::binary_and(a, b); };
  EnumerationNetwork::binary_t binary_or_func = [](const auto& a, const auto& b){ return kitty::binary_or(a, b); };
  EnumerationNetwork::ternary_t ternary_maj_func = [](const auto& a, const auto& b, const auto& c){ return kitty::ternary_majority(a, b, c); };
  operations.push_back(binary_and_func);

  create_and_try_functions( classes, xs, K, operations );

  /*printf("\n");
  printf("NOT CROSSING FREE\n");
  for( auto const& tt : classes )
  {
    kitty::print_binary( tt );
    printf("\n");
  }*/

  return 0;
}

#endif // MOCKTURTLE_PLANAR_FOUR_INPUT_HPP
