/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_RECIPROCAL_MATRIX_IMPL_H
#define OILAB_RECIPROCAL_MATRIX_IMPL_H

#include "RationalMatrix.h"
#include "../IO/Logger.h"
#include "../Math/BestRationalApproximation.h"
#include "../Math/IntegerMath.h"
#include <cfloat>  // FLT_EPSILON
#include <iomanip>

namespace oILAB {

/**********************************************************************/
template<int dim>
std::pair<typename RationalMatrix<dim>::MatrixDimI, typename RationalMatrix<dim>::IntScalarType> RationalMatrix<dim>::compute(
    const RationalMatrix<dim>::MatrixDimD& R)
{
    // Find the BestRationalApproximation of each entry
    MatrixDimI nums(MatrixDimI::Zero());
    MatrixDimI dens(MatrixDimI::Ones());

    IntScalarType sigma = 1;
    for(int i = 0; i < dim; ++i) {
        for(int j = 0; j < dim; ++j) {
            BestRationalApproximation bra(R(i, j), maxDen);
            nums(i, j) = bra.num;
            dens(i, j) = bra.den;
            sigma = IntegerMath<IntScalarType>::lcm(sigma, bra.den);
        }
    }

    MatrixDimI im(MatrixDimI::Zero());
    for(int i = 0; i < dim; ++i) {
        for(int j = 0; j < dim; ++j) {
            im(i, j) = nums(i, j) * (sigma / dens(i, j));
        }
    }

    const double error = (im.template cast<double>() / sigma - R).norm() / (dim * dim);
    if(error > FLT_EPSILON) {
        Logger::debug() << "maxDen=" << maxDen;
        Logger::debug() << "im=\n" << std::setprecision(15) << std::scientific << im.template cast<double>() / sigma;
        Logger::debug() << "= 1/" << sigma << "*\n" << std::setprecision(15) << std::scientific << im;
        Logger::debug() << "R=\n" << std::setprecision(15) << std::scientific << R;
        throw std::runtime_error("Rational Matrix failed, check maxDen");
    }

    return std::make_pair(im, sigma);
}

template<int dim>
std::pair<typename RationalMatrix<dim>::MatrixDimI, typename RationalMatrix<dim>::IntScalarType> RationalMatrix<dim>::reduce(
    const MatrixDimI& Rn, const MatrixDimI& Rd)
{
    if(Rd.any() == 0) throw std::runtime_error("Rational Matrix construction failed: denominator matrix has zeros");
    MatrixDimI im(MatrixDimI::Zero());
    MatrixDimI RnReduced(Rn);
    MatrixDimI RdReduced(Rd);

    // reduce the ratios
    for(int i = 0; i < dim; ++i) {
        for(int j = 0; j < dim; ++j) {
            RnReduced(i, j) = Rn(i, j) / IntegerMath<IntScalarType>::gcd(Rn(i, j), Rd(i, j));
            RdReduced(i, j) = Rd(i, j) / IntegerMath<IntScalarType>::gcd(Rn(i, j), Rd(i, j));
        }
    }
    IntScalarType sigma = IntegerMath<IntScalarType>::lcm(RdReduced.cwiseAbs());
    for(int i = 0; i < dim; ++i) {
        for(int j = 0; j < dim; ++j) {
            im(i, j) = RnReduced(i, j) * sigma / RdReduced(i, j);
        }
    }
    if(IntegerMath<IntScalarType>::gcd(IntegerMath<IntScalarType>::gcd(im.cwiseAbs()), sigma) != 1) {
        Logger::debug() << Rn;
        Logger::debug() << Rd;
        Logger::debug() << RnReduced;
        Logger::debug() << RdReduced;
        Logger::debug() << im;
        Logger::debug() << sigma;
        throw std::runtime_error("RationalMatrix<dim>::reduce failed: gcd!=1");
    }
    return std::make_pair(im, sigma);
}
/**********************************************************************/
template<int dim>
RationalMatrix<dim>::RationalMatrix(const MatrixDimD& R)
    : /* init */ returnPair(compute(R)),
      /* init */ integerMatrix(returnPair.first),
      /* init */ mu(returnPair.second)
{
}

template<int dim>
RationalMatrix<dim>::RationalMatrix(const MatrixDimI& Rn, const MatrixDimI& Rd)
try
    : /* init */ returnPair(reduce(Rn, Rd)),
      /* init */ integerMatrix(returnPair.first),
      /* init */ mu(returnPair.second) {
}
catch(std::runtime_error& e) {
    Logger::error() << e.what();
    throw(std::runtime_error("Rational Matrix construction failed. "));
}
template<int dim>
RationalMatrix<dim>::RationalMatrix(const MatrixDimI& Rn, const IntScalarType& Rd)
    : /* init */ returnPair(std::make_pair(Rn, Rd)),
      /* init */ integerMatrix(returnPair.first),
      /* init */ mu(returnPair.second)
{
}

template<int dim>
typename RationalMatrix<dim>::MatrixDimD RationalMatrix<dim>::asMatrix() const
{
    return integerMatrix.template cast<double>() / mu;
}

}  // namespace oILAB

#endif