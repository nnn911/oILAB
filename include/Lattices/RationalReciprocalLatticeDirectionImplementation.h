/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_RATIONAL_RECIPROCAL_LATTICE_DIRECTION_IMPL_H
#define OILAB_RATIONAL_RECIPROCAL_LATTICE_DIRECTION_IMPL_H

#include "RationalReciprocalLatticeDirection.h"

namespace oILAB {
/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim>::RationalReciprocalLatticeDirection(const Rational<IntScalarType>& _rat,
                                                                            const ReciprocalLatticeDirection<dim>& _dir)
    : /* init */ rat(_rat)
      /* init */,
      dir(_dir)
{
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim>::RationalReciprocalLatticeDirection(const Rational<IntScalarType>& _rat,
                                                                            const ReciprocalLatticeVector<dim>& v)
    : /* init */ RationalReciprocalLatticeDirection(_rat, ReciprocalLatticeDirection<dim>(v))
{
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim>::RationalReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v)
    : /* init */ rat(Rational<IntScalarType>(IntegerMath<IntScalarType>::gcd(v), 1)),
      /* init */ dir(v)
{
}

/**********************************************************************/

template<int dim>
typename RationalReciprocalLatticeDirection<dim>::VectorDimD RationalReciprocalLatticeDirection<dim>::cartesian() const
{
    return dir.cartesian() * rat.asDouble();
}

/**********************************************************************/
template<int dim>
Rational<typename RationalReciprocalLatticeDirection<dim>::IntScalarType> RationalReciprocalLatticeDirection<dim>::dot(
    const LatticeVector<dim>& other) const
{
    return rat * other.dot(dir);
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator*(const IntScalarType& scalar) const
{
    return RationalReciprocalLatticeDirection<dim>(rat * scalar, dir);
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator/(const IntScalarType& scalar) const
{
    return RationalReciprocalLatticeDirection<dim>(rat / scalar, dir);
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator+(
    const RationalReciprocalLatticeDirection<dim>& other) const
{
    assert(&dir.lattice == &other.dir.lattice && "Rational Lattice Vector Type belong to different Lattices.");
    const VectorDimI temp(rat.n * other.rat.d * dir.reciprocalLatticeVector() + other.rat.n * rat.d * other.dir.reciprocalLatticeVector());
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
    const ReciprocalLatticeVector<dim> v((temp / gcd).eval(), dir.lattice);
    return RationalReciprocalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), ReciprocalLatticeDirection<dim>(v));
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator-(
    const RationalReciprocalLatticeDirection<dim>& other) const
{
    assert(&dir.lattice == &other.dir.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
    const VectorDimI temp(rat.n * other.rat.d * dir.reciprocalLatticeVector() - other.rat.n * rat.d * other.dir.reciprocalLatticeVector());
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(temp));
    const ReciprocalLatticeVector<dim> v((temp / gcd).eval(), dir.lattice);
    return RationalReciprocalLatticeDirection<dim>(Rational(gcd, rat.d * other.rat.d), ReciprocalLatticeDirection<dim>(v));
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator+(const ReciprocalLatticeVector<dim>& other) const
{
    assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
    return this->operator+(
        RationalReciprocalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), ReciprocalLatticeDirection<dim>(other)));
}

/**********************************************************************/
template<int dim>
RationalReciprocalLatticeDirection<dim> RationalReciprocalLatticeDirection<dim>::operator-(const ReciprocalLatticeVector<dim>& other) const
{
    assert(&dir.lattice == &other.lattice && "ReciprocalLatticeVectorType belong to different Lattices.");
    const IntScalarType gcd(IntegerMath<IntScalarType>::gcd(other));
    return this->operator-(
        RationalReciprocalLatticeDirection<dim>(Rational<IntScalarType>(gcd, 1), ReciprocalLatticeDirection<dim>(other)));
}

/**********************************************************************/
template<int dim>
double RationalReciprocalLatticeDirection<dim>::squaredNorm() const
{
    return dir.reciprocalLatticeVector().squaredNorm() * std::pow(rat.asDouble(), 2);
}

template<int dim>
RationalReciprocalLatticeDirection<dim> operator*(const typename RationalReciprocalLatticeDirection<dim>::IntScalarType& scalar,
                                                  const RationalReciprocalLatticeDirection<dim>& L)
{
    return L * scalar;
}

}  // namespace oILAB

#endif