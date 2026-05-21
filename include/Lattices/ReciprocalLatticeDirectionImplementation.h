/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_RECIPROCAL_LATTICE_DIRECTION_IMPL_H
#define OILAB_RECIPROCAL_LATTICE_DIRECTION_IMPL_H

#include "../Math/RLLL.h"
#include "ReciprocalLatticeDirection.h"
#include <ostream>

namespace oILAB {
template<int dim>
ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v)
    : /* init */ ReciprocalLatticeVector<dim>(((v.squaredNorm() == 0) ? v : (v / abs(IntegerMath<IntScalarType>::gcd(v))).eval()),
                                              *v.lattice)
{
}

template<int dim>
ReciprocalLatticeDirection<dim>::ReciprocalLatticeDirection(const VectorDimI& v, const Lattice<dim>& lat)
    : /* base init */ ReciprocalLatticeVector<dim>(((v.squaredNorm() == 0) ? v : (v / IntegerMath<IntScalarType>::gcd(v)).eval()), lat)
{
}

template<int dim>
std::basic_ostream<char>& operator<<(std::basic_ostream<char>& s, const ReciprocalLatticeDirection<dim>& m)
{
    return s << m.reciprocalLatticeVector().transpose();
}

template<int dim>
double ReciprocalLatticeDirection<dim>::planeSpacing() const
{
    return 1.0 / cartesian().norm();
}

template<int dim>
int ReciprocalLatticeDirection<dim>::stacking() const
{
    RLLL rlll(this->lattice->latticeBasis, 0.75);
    auto structureMatrix = rlll.reducedBasis();
    auto U = rlll.unimodularMatrix();

    Lattice<dim> reducedLattice(structureMatrix);
    VectorDimI temp = U.transpose() * (*this);
    ReciprocalLatticeVector<dim> r(temp, reducedLattice);
    LatticeDirection<dim> vector(reducedLattice);
    vector = reducedLattice.latticeDirection(r.cartesian());
    return abs(vector.dot(r));
}

}  // namespace oILAB
#endif