/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_RATIONAL_LATTICE_DIRECTION_IMPL_H
#define OILAB_RATIONAL_LATTICE_DIRECTION_IMPL_H

#include "RationalLatticeDirection.h"

namespace oILAB {
/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim>::RationalLatticeDirection(const Rational<IntScalarType>& _rat, const LatticeDirection<dim>& _dir)
    : /* init */ rat(_rat)
      /* init */,
      dir(_dir)
{
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim>::RationalLatticeDirection(const Rational<IntScalarType>& _rat, const LatticeVector<dim>& v)
    : /* init */ RationalLatticeDirection(_rat, LatticeDirection<dim>(v))
{
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim>::RationalLatticeDirection(const LatticeVector<dim>& v)
    : /* init */ rat(Rational<IntScalarType>(IntegerMath<IntScalarType>::gcd(v), 1)),
      /* init */ dir(v)
{
}

/**********************************************************************/

template<int dim>
typename RationalLatticeDirection<dim>::VectorDimD RationalLatticeDirection<dim>::cartesian() const
{
    return dir.cartesian() * rat.asDouble();
}

/**********************************************************************/
template<int dim>
Rational<typename RationalLatticeDirection<dim>::IntScalarType> RationalLatticeDirection<dim>::dot(
    const ReciprocalLatticeVector<dim>& other) const
{
    return rat * dir.dot(other);
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator*(const IntScalarType& scalar) const
{
    return RationalLatticeDirection<dim>(rat * scalar, dir);
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator/(const IntScalarType& scalar) const
{
    return RationalLatticeDirection<dim>(rat / scalar, dir);
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator+(const RationalLatticeDirection<dim>& other) const
{
    assert(dir.lattice == other.dir.lattice && "Rational Lattice Vector Type belong to different Lattices.");
    // const VectorDimI temp(rat.n * other.rat.d * dir + other.rat.n * rat.d * other.dir);
    const VectorDimI temp(rat.n * other.rat.d * dir.latticeVector() + other.rat.n * rat.d * other.dir.latticeVector());
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
    const LatticeVector<dim> v((temp / gcd).eval(), *dir.lattice);
    return RationalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), LatticeDirection<dim>(v));
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator-(const RationalLatticeDirection<dim>& other) const
{
    assert(dir.lattice == other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
    // const VectorDimI temp(rat.n * other.rat.d * dir - other.rat.n * rat.d * other.dir);
    const VectorDimI temp(rat.n * other.rat.d * dir.latticeVector() - other.rat.n * rat.d * other.dir.latticeVector());
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
    const LatticeVector<dim> v((temp / gcd).eval(), *dir.lattice);
    return RationalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), LatticeDirection<dim>(v));
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator+(const LatticeVector<dim>& other) const
{
    assert(dir.lattice == other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
    return this->operator+(RationalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), LatticeDirection<dim>(other)));
}

/**********************************************************************/
template<int dim>
RationalLatticeDirection<dim> RationalLatticeDirection<dim>::operator-(const LatticeVector<dim>& other) const
{
    assert(dir.lattice == other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
    return this->operator-(RationalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), LatticeDirection<dim>(other)));
}

/**********************************************************************/
template<int dim>
double RationalLatticeDirection<dim>::squaredNorm() const
{
    // return dir.squaredNorm() * std::pow(rat.asDouble(), 2);
    return dir.latticeVector().squaredNorm() * std::pow(rat.asDouble(), 2);
}

template<int dim>
RationalLatticeDirection<dim> operator*(const typename RationalLatticeDirection<dim>::IntScalarType& scalar,
                                        const RationalLatticeDirection<dim>& L)
{
    return L * scalar;
}

}  // namespace oILAB

#endif