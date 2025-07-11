#pragma once
#ifndef _SYNCMER_H_
#define _SYNCMER_H_

#include "circularqueue.h"
#include "config.h"
#include "lcptools/lcptools_ho.h"
#include <concepts>
#include <gtl/bit_vector.hpp>
#include <limits>
#include <stdint.h>
#include <type_traits>

namespace syncmer {

template <class INT_TYPE> struct murmur {
  using type = INT_TYPE;
};
template <> struct murmur<uint64_t> {
  using type = uint64_t;
  type operator()(type val) {
    type v1 = MurmurHash3_32(&val, sizeof(val), 1);
    v1 ^= +0x9e3779b9 + (v1 << 6) + (v1 >> 2);
    return (v1 << 32) | MurmurHash3_32(&val, sizeof(val), v1);
  }
};

template <> struct murmur<uint32_t> {
  using type = uint32_t;
  type operator()(type val) {
    type v1 = MurmurHash3_32(&val, sizeof(val), 1);
    return v1;
  }
};
template <> struct murmur<uint16_t> {
  using type = uint32_t;
  type operator()(type val) {
    type v1 = MurmurHash3_32(&val, sizeof(val), 1);
    return v1;
  }
};
template <class Rule>
concept SyncmerRule =
    std::is_integral_v<typename Rule::int_type> && requires(Rule k) {
      { typename Rule::int_type{} } -> std::same_as<typename Rule::int_type>;
      { Rule::K } -> std::convertible_to<int32_t>;
      { Rule::S } -> std::convertible_to<int32_t>;
      { Rule::T } -> std::convertible_to<int32_t>;
      {
        Rule::hash(typename Rule::int_type{})
      } -> std::convertible_to<typename Rule::int_type>;
      {
        Rule::is_syncmer(typename Rule::int_type{})
      } -> std::convertible_to<bool>;
      {
        Rule::is_syncmer(smq<typename Rule::int_type, Rule::K>{})
      } -> std::convertible_to<bool>;
    };

#define CLOSED_T -1
template <int _K, int _S, int _T, typename INT_TYPE = uint32_t,
          class Hash = murmur<INT_TYPE>>
struct OpenSyncmerRule {
  using int_type = INT_TYPE;

  static constexpr int K = _K;
  static constexpr int S = _S;
  static constexpr int T = _T;
   static constexpr int_type K_MASK = (1 << K) - 1;
  static constexpr int_type S_MASK = (1 << S) - 1;
    static constexpr int_type hash(int_type a) noexcept { return Hash{}(a); }
  static constexpr bool is_syncmer(int_type a) noexcept {
    int min_index = -1;
    int_type min_val = std::numeric_limits<int_type>::max();
    for (int i = 0; i < K - S; ++i) {
      int_type masked_s = (a >> i) & S_MASK;
      int_type hashed_s = hash(masked_s);
      if (hashed_s < min_val) {
        min_index = i;
        min_val = hashed_s;
      }
    }
    return T == min_index;
  }
  static constexpr bool is_syncmer(const smq<int_type, K> &q) noexcept {
    return q.min() == q[T];
  }
};

template <int _K, int _S, int _TIGNORED, typename INT_TYPE = uint32_t,
          class Hash = murmur<INT_TYPE>>
struct ClosedSyncmerRule {
  using int_type = INT_TYPE;
  static constexpr int K = _K;
  static constexpr int S = _S;
  static constexpr int T = CLOSED_T;
  static constexpr int_type hash(int_type a) noexcept { return Hash{}(a); }
  static constexpr int_type K_MASK = (1 << K) - 1;
  static constexpr int_type S_MASK = (1 << S) - 1;
  static constexpr bool is_syncmer(uint64_t a) noexcept {
    int min_index = -1;
    int_type min_val = std::numeric_limits<int_type>::max();
    for (int i = 0; i < K - S; ++i) {
      int_type masked_s = (a >> i) & S_MASK;
      int_type hashed_s = hash(masked_s);
      if (hashed_s < min_val) {
        min_index = i;
        min_val = hashed_s;
      }
    }
    return min_index == 0 || min_index == K - S;
  }
  static constexpr bool is_syncmer(const smq<int_type, K> &q) noexcept {
    return q.min() == q.front() || q.min() == q.back();
  }
};

template <int _K, int Numerator, int Denominator, typename INT_TYPE = uint32_t,
          class Hash = murmur<INT_TYPE>>
struct DensitySyncmerRule {
  using int_type = INT_TYPE;
  static constexpr int K = _K;
  static constexpr int S = Numerator;
  static constexpr int T = Denominator;
  static constexpr double density =
      (std::numeric_limits<int_type>::max() / double(Denominator)) * Numerator;
  static int_type hash(int_type a) noexcept { return Hash{}(a); }
  static bool is_syncmer(int_type a) noexcept { return hash(a) < density; }
  static bool is_syncmer(const smq<int_type, _K> &) noexcept { return false; }
};

template <SyncmerRule Rule> struct syncmer_precomputed { // TODO maybe bloom filter version
  using rule = Rule;
  const size_t size = 1 << (2 * Rule::K - 1);
  gtl::bit_vector lut;
  syncmer_precomputed() : lut(size) {
    for (size_t i = 0; i < size; ++i) {
      lut.set(i, rule::is_syncmer(i & KMER_SPACE_MASK));
    }
  }
  constexpr bool operator[](size_t idx) const noexcept { return lut[idx]; }
};

} // namespace syncmer


#endif 