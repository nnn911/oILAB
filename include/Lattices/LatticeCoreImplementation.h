/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_LATTICE_CORE_IMPL_H
#define OILAB_LATTICE_CORE_IMPL_H

#include "../../include/Lattices/LatticeCore.h"
#include "../../include/IO/Logger.h"
#include "../../include/Math/BestRationalApproximation.h"
#include <Eigen/Dense>

namespace oILAB {

template<int dim>
typename LatticeCore<dim>::VectorDimI LatticeCore<dim>::rationalApproximation(VectorDimD nd)
{
    Eigen::Array<IntScalarType, dim, 1> nums = VectorDimI::Zero();

    if(nd.squaredNorm() > 0.0) {
        const Eigen::Array<double, dim, 1> nda(
            nd.array().abs());  // vector of close-to-integer numbers corresponding to lattice coordinates
        size_t maxID = 0;
        const double maxVal(nda.maxCoeff(&maxID));
        nd /= maxVal;  // make each value of nd in [-1:1]

        nums = VectorDimI::Ones();
        Eigen::Array<IntScalarType, dim, 1> dens = VectorDimI::Ones();
        IntScalarType denProd = 1;

        for(int k = 0; k < dim; ++k) {
            // BestRationalApproximation bra(nd(k), 10000);
            BestRationalApproximation bra(nd(k), 1000);

            nums(k) = bra.num;
            dens(k) = bra.den;
            denProd *= bra.den;
        }

        for(int k = 0; k < dim; ++k) {
            nums(k) *= (denProd / dens(k));
        }
    }

    return nums.matrix();
}

template<int dim>
typename LatticeCore<dim>::VectorDimI LatticeCore<dim>::integerCoordinates(const VectorDimD& d, const MatrixDimD& invA)
{
    const VectorDimD nd(invA * d);
    const VectorDimD rd(nd.array().round());
    if((nd - rd).norm() > roundTol) {
        Logger::debug() << "nd=" << nd.transpose();
        Logger::debug() << "rd=" << rd.transpose();
        Logger::debug() << "rounding error = |nd-rd| = " << (nd - rd).norm();
        throw(std::runtime_error("Input vector is not a lattice vector"));
    }
    return rd.template cast<IntScalarType>();
}

}  // namespace oILAB

#endif