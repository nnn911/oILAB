/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_LATTICE_VECTOR_IMPL_H
#define OILAB_LATTICE_VECTOR_IMPL_H

#include "LatticeVector.h"

namespace oILAB {

/**********************************************************************/
template<int dim>
typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base()
{
    return *this;
}

/**********************************************************************/
template<int dim>
const typename LatticeVector<dim>::BaseType& LatticeVector<dim>::base() const
{
    return *this;
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>::LatticeVector()
    : /* init */ BaseType(VectorDimI::Zero()),
      /* init */ lattice(nullptr)
{
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>::LatticeVector(const Lattice<dim>& lat)
    : /* init */ BaseType(VectorDimI::Zero()),
      /* init */ lattice(&lat)
{
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>::LatticeVector(const VectorDimD& d, const Lattice<dim>& lat)
    : /* init */ BaseType(LatticeCore<dim>::integerCoordinates(d, lat.reciprocalBasis.transpose())),
      /* init */ lattice(&lat)
{ /*!@param[in] d vector in real space
   * Constructs *this by mapping d to the lattice
   */
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>::LatticeVector(const VectorDimI& other, const Lattice<dim>& lat)
    : /* init base */ BaseType(other),
      /* init      */ lattice(&lat)
{
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>& LatticeVector<dim>::operator=(const LatticeVector<dim>& other)
{
    if(lattice == nullptr) lattice = other.lattice;
    else assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    base() = other.base();
    return *this;
}
/**********************************************************************/
template<int dim>
LatticeVector<dim>& LatticeVector<dim>::operator=(LatticeVector<dim>&& other)
{
    if(lattice == nullptr) lattice = other.lattice;
    else assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    base() = other.base();
    return *this;
}

/**********************************************************************/
template<int dim>
LatticeVector<dim> LatticeVector<dim>::operator+(const LatticeVector<dim>& other) const
{
    assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    VectorDimI temp = static_cast<VectorDimI>(*this) + static_cast<VectorDimI>(other);
    return LatticeVector<dim>(temp, *lattice);
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>& LatticeVector<dim>::operator+=(const LatticeVector<dim>& other)
{
    assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    base() += other.base();
    return *this;
}

/**********************************************************************/
template<int dim>
LatticeVector<dim> LatticeVector<dim>::operator-(const LatticeVector<dim>& other) const
{
    assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    VectorDimI temp = static_cast<VectorDimI>(*this) - static_cast<VectorDimI>(other);
    return LatticeVector<dim>(temp, *lattice);
}

/**********************************************************************/
template<int dim>
LatticeVector<dim>& LatticeVector<dim>::operator-=(const LatticeVector<dim>& other)
{
    assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    base() -= other.base();
    return *this;
}

/**********************************************************************/
template<int dim>
LatticeVector<dim> LatticeVector<dim>::operator*(const LatticeVector<dim>::IntScalarType& scalar) const
{
    VectorDimI temp = static_cast<VectorDimI>(*this) * scalar;
    return LatticeVector<dim>(temp, *lattice);
}
/**********************************************************************/
template<int dim>
typename LatticeVector<dim>::IntScalarType LatticeVector<dim>::dot(const ReciprocalLatticeVector<dim>& other) const
{
    assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
}

/**********************************************************************/
template<int dim>
typename LatticeVector<dim>::IntScalarType LatticeVector<dim>::dot(const ReciprocalLatticeDirection<dim>& other) const
{
    assert(lattice == other.lattice && "LatticeVectors belong to different Lattices.");
    return dot(other.reciprocalLatticeVector());
}

/**********************************************************************/
template<int dim>
typename LatticeVector<dim>::VectorDimD LatticeVector<dim>::cartesian() const
{
    return lattice->latticeBasis * this->template cast<double>();
}

/**********************************************************************/
template<int dim>
LatticeVector<dim> operator*(const typename LatticeVector<dim>::IntScalarType& scalar, const LatticeVector<dim>& L)
{
    return L * scalar;
}

template<int dim>
LatticeVector<dim> operator*(const int& scalar, const LatticeVector<dim>& L)
{
    return L * scalar;
}

}  // namespace oILAB

#endif