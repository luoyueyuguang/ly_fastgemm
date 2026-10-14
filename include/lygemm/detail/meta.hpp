// lygemm -- compile time helpers.
// SPDX-License-Identifier: MIT
#ifndef LYGEMM_META_HPP
#define LYGEMM_META_HPP

#include <cstddef>
#include <type_traits>
#include <utility>

#include "lygemm/config.hpp"

namespace lygemm {
namespace meta {

template <std::size_t A, std::size_t B>
struct min : std::integral_constant<std::size_t, (A < B ? A : B)> {};

template <std::size_t A, std::size_t B>
struct max : std::integral_constant<std::size_t, (A > B ? A : B)> {};

template <std::size_t V, std::size_t Lo, std::size_t Hi>
struct clamp
    : std::integral_constant<std::size_t, (V < Lo ? Lo : (V > Hi ? Hi : V))> {};

/* Calls f(std::integral_constant<std::size_t, I>{}) for every I in [0, N).
 * The fold expression makes the calls compile time constants, which is what
 * lets the kernels above be fully unrolled and the tile sizes be template
 * arguments of the micro kernels.
 *
 * Both the fold and the calls it makes have to be inlined: if the compiler is
 * allowed to outline an unrolled body, the fragment it accumulates into is
 * passed by reference and its register tile is demoted to the stack, which
 * costs several times the whole kernel. */
template <class F, std::size_t... I>
LYGEMM_ALWAYS_INLINE void static_for_impl(F&& f, std::index_sequence<I...>) {
  (f(std::integral_constant<std::size_t, I>{}), ...);
}

template <std::size_t N, class F>
LYGEMM_ALWAYS_INLINE void static_for(F&& f) {
  static_for_impl(static_cast<F&&>(f), std::make_index_sequence<N>{});
}

} // namespace meta
} // namespace lygemm

#endif /* LYGEMM_META_HPP */
