// -*- C++ -*-
#pragma once

#include "vector.h"

template<int N, typename T>
using tensor_t = vector_t<N, vector_t<N, T> >;
