// -*- C++ -*-
#pragma once

#include <array>
#include <iterator>
#include <ostream>
#include <cmath>

template<int N, typename T>
struct vector_t {
  template<typename... Args> vector_t(Args const &... args) : data{args...} {}
  template<typename I> T const &operator[](I const &i) const { return data[i]; }
  template<typename I> T &operator[](I const &i) { return data[i]; }
  std::array<T, N> data;
};

template<int N, typename T>
std::ostream &operator<<(std::ostream &output, vector_t<N, T> const &obj) {
  if (N > 0) {
    output << "[" << obj[0];
    for (auto i = 1; i != N; ++i) output << ", " << obj[i];
    output << "]";
  }
  return output;
}

template<int N, typename T>
auto operator+(vector_t<N, T> const &obj, vector_t<N, T> const &rhs) {
  vector_t<N, T> lhs;
  for (auto i = 0; i != N; ++i) lhs[i] = obj[i] + rhs[i];
  return lhs;
}

template<int N, typename T>
auto operator-(vector_t<N, T> const &obj, vector_t<N, T> const &rhs) {
  vector_t<N, T> lhs;
  for (auto i = 0; i != N; ++i) lhs[i] = obj[i] - rhs[i];
  return lhs;
}

template<int N, typename T>
auto operator+(vector_t<N, T> const &obj, T const &rhs) {
  vector_t<N, T> lhs;
  for (auto i = 0; i != N; ++i) lhs[i] = obj[i] + T(rhs);
  return lhs;
}

template<int N, typename T>
auto operator-(vector_t<N, T> const &obj, T const &rhs) {
  vector_t<N, T> lhs;
  for (auto i = 0; i != N; ++i) lhs[i] = obj[i] - T(rhs);
  return lhs;
}

template<int N, typename T>
auto operator/(vector_t<N, T> const &obj, T const &rhs) {
  vector_t<N, T> lhs;
  for (auto i = 0; i != N; ++i) lhs[i] = obj[i] / rhs;
  return lhs;
}

template<int N, typename T>
auto operator*(vector_t<N, T> const &obj, T const &rhs) {
  vector_t<N, T> lhs;
  for (auto i = 0; i != N; ++i) lhs[i] = obj[i] * rhs;
  return lhs;
}

template<int N, typename T, typename RHS>
auto dot(vector_t<N, T> const &obj, RHS &&rhs) {
  T lhs{0};
  for (auto i = 0; i != N; ++i) lhs += obj[i] * rhs[i];
  return lhs;
}

template<int N, typename T>
auto cross(vector_t<N, T> const &obj, vector_t<N, T> const &rhs) {
  vector_t<N, T> lhs{obj[1]*rhs[2]-obj[2]*rhs[1],
                     obj[2]*rhs[0]-obj[0]*rhs[2],
                     obj[0]*rhs[1]-obj[1]*rhs[0]};
  return lhs;
}

template<int N, typename T>
auto mag2(vector_t<N, T> const &obj) { return dot(obj, obj); }

template<int N, typename T>
auto mag(vector_t<N, T> const &obj) { return std::sqrt(mag2(obj)); }
