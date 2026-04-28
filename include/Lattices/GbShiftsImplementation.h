/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_GB_SHIFTS_IMPL_H
#define OILAB_GB_SHIFTS_IMPL_H

#include "GbShifts.h"
#include "../Utilities/randomInteger.h"
#include "../IO/Logger.h"

namespace oILAB {

template<int dim>
GbShifts<dim>::GbShifts(const Gb<dim>& gb,
                        const ReciprocalLatticeVector<dim>& axis,
                        const std::vector<LatticeVector<dim>>& gbCslVectors,
                        const double& bhalfMax)
    : gb(gb), axis(axis), gbCslVectors(gbCslVectors), bShiftPairs(getbShiftPairs(gb, gbCslVectors, bhalfMax))
{
    Logger::info() << "--------------------GBShifts class construction---------------------------";
    Logger::info() << "GB CSL vectors =";
    for(const auto& elem : gbCslVectors) Logger::info() << elem.cartesian().transpose();

    Logger::info() << "GB reciprocal CSL vectors =";
    Eigen::Matrix<double, dim, dim - 1> gbCslBasis;
    for(int i = 0; i < dim - 1; ++i) gbCslBasis.col(i) = gbCslVectors[i].cartesian();
    Eigen::Matrix<double, dim, dim - 1> gbCslReciprocalBasis = gbCslBasis.completeOrthogonalDecomposition().pseudoInverse().transpose();
    Logger::info() << gbCslReciprocalBasis.transpose();

    Logger::info() << "Maximum b < " << 2 * bhalfMax * gb.bc.A.latticeBasis.col(0).norm();

    VectorDimD normal;
    if(dim == 3)
        normal = gbCslVectors[0].cross(gbCslVectors[1]).cartesian().normalized();
    else
        normal = gbCslVectors[0].cross().cartesian().normalized();

    Logger::info() << "Exploring the following translation-shift pairs:";
    for(const auto& [b, s] : bShiftPairs) {
        Logger::info() << "b = " << b.cartesian().transpose() << "; s = " << s.transpose();
        // if (abs(s.dot(normal)) > FLT_EPSILON)
        if(abs(s.dot(normal)) > 1e-6) throw std::runtime_error("GBShifts construction failed - shifts are not parallel to the GB.");

        Eigen::MatrixXd shiftCoordinates = gbCslReciprocalBasis.transpose() * s;
        if((shiftCoordinates.array() < -FLT_EPSILON).any() || (shiftCoordinates.array() > 1 + FLT_EPSILON).any()) {
            Logger::debug() << "Shift coordinates = " << shiftCoordinates.transpose();
            throw std::runtime_error("GB shifts are not in the area spanned by the GB CSL vectors.");
        }
    }
    Logger::info() << "----------------------------";
}

template<int dim>
std::vector<std::pair<LatticeVector<dim>, typename GbShifts<dim>::VectorDimD>> GbShifts<dim>::getbShiftPairs(
    const Gb<dim>& gb, const std::vector<LatticeVector<dim>>& gbCslVectors, const double& bhalfMax)
{
    std::vector<std::pair<LatticeVector<dim>, VectorDimD>> output;

    assert(gbCslVectors.size() == dim - 1);
    auto nC = gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());

    // form a CSL cell for modulo operations
    auto gbPlaneParallelCslBasis = gb.bc.csl.planeParallelLatticeBasis(nC, true);
    std::vector<LatticeVector<dim>> cslSubLatticeVectors;
    cslSubLatticeVectors.push_back(gbPlaneParallelCslBasis[0].latticeVector());
    cslSubLatticeVectors.push_back(gbCslVectors[0]);
    cslSubLatticeVectors.push_back(gbCslVectors[1]);
    auto cslPoints = gb.bc.csl.box(cslSubLatticeVectors);

    // form a T lattice cell for exploring translations
    auto nT = gb.getReciprocalLatticeDirectionInT(
        gb.bc.getReciprocalLatticeDirectionInD(gb.nA.reciprocalLatticeVector()).reciprocalLatticeVector());
    auto planeParallelBasisT = gb.T.planeParallelLatticeBasis(nT, true);
    std::vector<LatticeVector<dim>> latticeVectorsT;

    double latticeConstant = gb.bc.A.latticeBasis.col(0).norm();
    for(int i = 0; i < dim; ++i) {
        // scale the lattice vectors of T based on the bhalfMax parameter
        // int factor= floor(bhalfMax*latticeConstant/planeParallelBasisT[i].latticeVector().cartesian().norm()+FLT_EPSILON);
        int factor = floor(35 * bhalfMax * latticeConstant / planeParallelBasisT[i].latticeVector().cartesian().norm() + FLT_EPSILON);
        factor = (factor > 0 ? factor : 1);
        latticeVectorsT.push_back(factor * planeParallelBasisT[i].latticeVector());
    }
    auto points = gb.T.box(latticeVectorsT, "T.txt");

    VectorDimD shiftT, shiftC;
    shiftT << -0.5, -0.5, -0.5;
    shiftC << -0.5, -FLT_EPSILON, -FLT_EPSILON;
    for(auto& point : points) {
        // if (point.cartesian().norm() > bhalfMax*gb.bc.A.latticeBasis.col(0).norm())
        //     continue;
        LatticeVector<dim>::modulo(point, latticeVectorsT, shiftT);
        if(point.cartesian().norm() > bhalfMax * gb.bc.A.latticeBasis.col(0).norm()) continue;
        auto cslShift = LatticeVector<dim>((gb.bc.LambdaA * gb.basisT * point).eval(), gb.bc.dscl);
        for(const auto& cslPoint : cslPoints) {
            // only for those CSL points on the boundary
            if(gb.bc.getLatticeVectorInA(cslPoint).dot(gb.nA) != 0) continue;
            VectorDimD cslShiftCentered = cslShift.cartesian() + cslPoint.cartesian() - point.cartesian() / 2;
            LatticeVector<dim>::modulo(cslShiftCentered, cslSubLatticeVectors, shiftC);
            output.push_back(std::make_pair(point, cslShiftCentered));
        }
    }
    return output;
}
/*
template<int dim>
std::vector<LatticeVector<dim>> GbShifts<dim>::getGbCslVectors(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis)
{
    std::vector<LatticeVector<dim>> output;
    LatticeVector<dim> axisA(gb.bc.A.latticeDirection(axis.cartesian()).latticeVector());
    LatticeVector<dim> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
    for(int i=0; i<dim-1; ++i) {
        if (i==0) output.push_back(gb.getPeriodVector(axis));
        if (i==1) output.push_back(axisC);
    }
    return output;
}
 */

}  // namespace oILAB

#endif