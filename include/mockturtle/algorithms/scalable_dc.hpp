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
 \file window_utils.hpp
 \brief Utilities to collect small-scale sets of nodes

 \author Benjamin Hien
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <optional>
#include <set>
#include <type_traits>
#include <vector>

namespace mockturtle
{

template<typename Ntk>
class create_dc_windows_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  struct window
  {
    std::vector<node> inputs;
    std::vector<node> nodes;
    std::vector<signal> outputs;
  };

protected:
  /* constant node used to denote invalid window elements */
  static constexpr node INVALID_NODE{ 0 };

  static constexpr uint32_t NUM_ITERATIONS{ 5 };

public:
  explicit create_dc_windows_impl( Ntk const& ntk )
      : ntk( ntk ), path( ntk.size() ), refs( ntk.size() )
  {
  }

  void resize( uint32_t size )
  {
    path.resize( size );
    refs.resize( size );
  }

  void run( std::vector<node> const& leaves, uint32_t num_levels )
  {
    windows.clear(); // clear any previous state

    for ( const auto& leaf : leaves )
    {
      window w;

      // Dummy window content
      w.inputs.push_back( leaf );
      w.nodes.push_back( leaf );

      // Create a constant false signal as dummy output
      w.outputs.push_back( ntk.get_constant( false ) );

      windows.push_back( std::move( w ) );
    }
  }

protected:
  Ntk const& ntk;
  std::vector<node> visited;
  std::vector<node> path;
  std::vector<uint32_t> refs;
  std::vector<std::vector<node>> levels;

  std::vector<window> windows;
};

} // namespace mockturtle