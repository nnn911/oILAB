/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_LATTICE_DIRECTION_IMPL_H
#define OILAB_LATTICE_DIRECTION_IMPL_H

#include "LatticeDirection.h"

namespace oILAB {
template<int dim>
LatticeDirection<dim>::LatticeDirection(const LatticeVector<dim>& v)
    : /* base init */ LatticeVector<dim>(((v.squaredNorm() == 0) ? v : (v / IntegerMath<IntScalarType>::gcd(v)).eval()), v.lattice)
{
}

template<int dim>
LatticeDirection<dim>::LatticeDirection(const VectorDimI& v, const Lattice<dim>& lat)
    : /* base init */ LatticeVector<dim>(((v.squaredNorm() == 0) ? v : (v / IntegerMath<IntScalarType>::gcd(v)).eval()), lat)
{
}

template<int dim>
basic_ostream<char>& operator<<(basic_ostream<char>& s, const LatticeDirection<dim>& m)
{
    return s << m.latticeVector().transpose();
}

}  // namespace oILAB
#endif