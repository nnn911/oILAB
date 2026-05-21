/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef OILAB_BICRYSTAL_IMPL_H
#define OILAB_BICRYSTAL_IMPL_H

#include "../IO/Logger.h"
#include "BiCrystal.h"
#include <numbers>
#include <fstream>
#include "../Math/RLLL.h"
#include "../Utilities/Rotation.h"
#include <bitset>
#include <optional>
#include <set>
#include <array>

namespace oILAB {

template<int dim>
typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getM(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd)
{
    typename BiCrystal<dim>::MatrixDimI M = BiCrystal<dim>::MatrixDimI::Identity();
    for(int i = 0; i < dim; ++i) {
        const auto& dii(sd.matrixD()(i, i));
        M(i, i) = dii / IntegerMath<IntScalarType>::gcd(rm.mu, dii);
    }
    return M;
}

template<int dim>
typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getN(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd)
{
    typename BiCrystal<dim>::MatrixDimI N = BiCrystal<dim>::MatrixDimI::Identity();
    for(int i = 0; i < dim; ++i) {
        const auto& dii = sd.matrixD()(i, i);
        N(i, i) = rm.mu / IntegerMath<IntScalarType>::gcd(rm.mu, dii);
    }
    return N;
}

template<int dim>
typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getLambdaA(const typename BiCrystal<dim>::MatrixDimI& M,
                                                               const typename BiCrystal<dim>::MatrixDimI& N)
{
    BiCrystal<dim>::MatrixDimI LambdaA;
    for(int col = 0; col < dim; ++col) {
        BiCrystal<dim>::VectorDimI x, y;
        for(int i = 0; i < dim; ++i) {
            IntScalarType a = N(i, i);
            IntScalarType b = -M(i, i);
            IntScalarType c = -(i == col);
            IntegerMath<IntScalarType>::solveDiophantine2vars(a, b, c, x(i), y(i));
        }
        LambdaA.col(col) = M * y;
    }
    return LambdaA;
}

template<int dim>
typename BiCrystal<dim>::MatrixDimI BiCrystal<dim>::getLambdaB(const typename BiCrystal<dim>::MatrixDimI& M,
                                                               const typename BiCrystal<dim>::MatrixDimI& N)
{
    typename BiCrystal<dim>::MatrixDimI LambdaB;
    for(int col = 0; col < dim; ++col) {
        BiCrystal<dim>::VectorDimI x, y;
        for(int i = 0; i < dim; ++i) {
            IntScalarType a = N(i, i);
            IntScalarType b = -M(i, i);
            IntScalarType c = (i == col);
            IntegerMath<IntScalarType>::solveDiophantine2vars(a, b, c, x(i), y(i));
        }
        LambdaB.col(col) = N * x;
    }
    return LambdaB;
}

template<int dim>
typename BiCrystal<dim>::MatrixDimD BiCrystal<dim>::getCSLBasis(const Lattice<dim>& A,
                                                                const Lattice<dim>& B,
                                                                const SmithDecomposition<dim>& sd,
                                                                const typename BiCrystal<dim>::MatrixDimI& M,
                                                                const typename BiCrystal<dim>::MatrixDimI& N)
{
    // The transition matrix is T=P/sigma, where P=rm.integerMatrix is
    // an integer matrix and sigma=rm.sigma is an integer
    // The integer matrix P can be decomposed as P=X*D*Y using the Smith decomposition.
    // X and Y are unimodular, and D is diagonal with D(k,k) dividing D(k+1,k+1)
    // The decomposition also computes the matices U and V such that D=U*P*V
    // From T=inv(A)*B=P/sigma=X*D*Y/sigma=X*D*inv(V)/sigma, we have
    // B1*(sigma*I)=A1*D
    // where
    // B1=B*V
    // A1=A*X
    // Since V and X are unimodular matrices, B1 and A1 are new bases
    // of the lattices B and A, respectively. Moreover, since
    // (sigma*I) and D are diagonal, the columns of B1 and A1 are
    // proportional, with rational proportionality factors different for each column.
    // For example, columns "i" read
    // b1_i*sigma=a1_i*D(i,i)
    // Therefore, the i-th primitive vectors of the CSL is
    // c_i=b1_i*sigma/gcd(sigma,D(i,i))=a1_i*D(i,i)/gcd(sigma,D(i,i))
    // or, in matrix form
    // C=B1*N=A1*M, that is
    // C=B*V*N=A*X*M
    // where M=diag(D(i,i)/gcd(sigma,D(i,i))) and
    //       N=diag(sigma/gcd(sigma,D(i,i))) and

    const auto C1(A.latticeBasis * (sd.matrixX() * M).template cast<double>());
    const auto C2(B.latticeBasis * (sd.matrixV() * N).template cast<double>());
    if((C1 - C2).norm() / C1.norm() > FLT_EPSILON || (C1 - C2).norm() / C2.norm() > FLT_EPSILON) {
        throw std::runtime_error("CSL calculation failed.");
    }

    /*
    if(useRLLL)
    {
        return RLLL(0.5*(C1+C2),0.75).reducedBasis();
    }
    else
    {
        return 0.5*(C1+C2);
    }
     */
    return 0.5 * (C1 + C2);
}

template<int dim>
typename BiCrystal<dim>::MatrixDimD BiCrystal<dim>::getDSCLBasis(const Lattice<dim>& A,
                                                                 const Lattice<dim>& B,
                                                                 const SmithDecomposition<dim>& sd,
                                                                 const typename BiCrystal<dim>::MatrixDimI& M,
                                                                 const typename BiCrystal<dim>::MatrixDimI& N)
{
    const auto D1(A.latticeBasis * sd.matrixX().template cast<double>() * N.template cast<double>().inverse());
    const auto D2(B.latticeBasis * sd.matrixV().template cast<double>() * M.template cast<double>().inverse());
    if((D1 - D2).norm() / D1.norm() > FLT_EPSILON || (D1 - D2).norm() / D2.norm() > FLT_EPSILON) {
        throw std::runtime_error("DSCL calculation failed.");
    }
    /*
    if(useRLLL)
    {
        return RLLL(0.5*(D1+D2),0.75).reducedBasis();
    }
    else
    {
        return 0.5*(D1+D2);
    }
     */
    return 0.5 * (D1 + D2);
}

template<int dim>
BiCrystal<dim>::BiCrystal(const Lattice<dim>& A_in, const Lattice<dim>& B_in, const bool& useRLLL)
try : /* init */ RationalMatrix
    <dim>(A_in.reciprocalBasis.transpose() * B_in.latticeBasis)
        /* init */,
        SmithDecomposition<dim>(this->integerMatrix)
        /* init */,
        A(A_in)
        /* init */,
        B(B_in)
        /* init */,
        M(getM(*this, *this))
        /* init */,
        N(getN(*this, *this))
        /* init */,
        sigmaA(round(M.template cast<double>().determinant()))
        /* init */,
        sigmaB(round(N.template cast<double>().determinant()))
        /* init */,
        sigma(std::abs(sigmaA) == std::abs(sigmaB) ? std::abs(sigmaA) : 0)
        /* init */,
        cslp(getCSLBasis(A, B, *this, M, N), MatrixDimD::Identity())
        /* init */,
        dsclp(getDSCLBasis(A, B, *this, M, N), MatrixDimD::Identity())
        /* init */,
        csl2cslp(useRLLL ? RLLL(cslp.latticeBasis, 0.75).unimodularMatrix() : MatrixDimI::Identity())
        /* init */,
        dscl2dsclp(useRLLL ? RLLL(dsclp.latticeBasis, 0.75).unimodularMatrix() : MatrixDimI::Identity())
        /* init */,
        csl(cslp.latticeBasis* csl2cslp.template cast<double>())
        /* init */,
        dscl(dsclp.latticeBasis* dscl2dsclp.template cast<double>())
        /* init */,
        Ap(A.latticeBasis* this->matrixX().template cast<double>())
        /* init */,
        Bp(B.latticeBasis* this->matrixV().template cast<double>())
        /* init */,
        LambdaA(getLambdaA(M, N))
        /* init */,
        LambdaB(getLambdaB(M, N))
    {
        if(true) {  // verify that CSL can be obtained as multiple of A and B

            const MatrixDimD tempA(A.reciprocalBasis.transpose() * csl.latticeBasis);
            if((tempA - tempA.array().round().matrix()).norm() / tempA.norm() > FLT_EPSILON) {
                // std::cout << (tempA-tempA.array().round().matrix()).norm()/tempA.norm() << std::endl;
                throw std::runtime_error("CSL is not a multiple of lattice A.");
            }

            const MatrixDimD tempB(B.reciprocalBasis.transpose() * csl.latticeBasis);
            if((tempB - tempB.array().round().matrix()).norm() / tempB.norm() > FLT_EPSILON) {
                // std::cout << (tempB-tempB.array().round().matrix()).norm()/tempB.norm() << std::endl;
                throw std::runtime_error("CSL is not a multiple of lattice B");
            }
        }

        if(true) {  // verify that A and B are multiples of DSCL

            const MatrixDimD tempA(dscl.reciprocalBasis.transpose() * A.latticeBasis);
            if((tempA - tempA.array().round().matrix()).norm() / tempA.norm() > FLT_EPSILON) {
                // std::cout << (tempA-tempA.array().round().matrix()).norm()/tempA.norm() << std::endl;
                throw std::runtime_error("Lattice A is not a multiple of the DSCL");
            }

            const MatrixDimD tempB(dscl.reciprocalBasis.transpose() * B.latticeBasis);
            if((tempB - tempB.array().round().matrix()).norm() / tempB.norm() > FLT_EPSILON) {
                // std::cout << (tempB-tempB.array().round().matrix()).norm()/tempB.norm() << std::endl;
                throw std::runtime_error("Lattice B is not a multiple of the DSCL");
            }
        }

        if(true) {  // verify LambdaA + LambdaB = I
            if(!(LambdaA + LambdaB).isIdentity()) {
                throw std::runtime_error("LambdaA + LambdaB != I");
            }
        }
    }

catch(std::runtime_error& e) {
    Logger::error() << e.what();
    throw(std::runtime_error("Bicrystal construction failed."));
}

template<int dim>
LatticeVector<dim> BiCrystal<dim>::getLatticeVectorInA(const LatticeVector<dim>& v) const
{
    VectorDimI integerCoordinates;

    if(v.lattice == &(this->A))
        return v;
    else if(v.lattice == &(this->csl))
        // U*M*csl2cslp*v
        integerCoordinates = this->matrixX() * M * csl2cslp * v;
    else
        throw(
            std::runtime_error("The input lattice vector should belong "
                               "to lattice A or the CSL"));
    auto temp = LatticeVector<dim>(integerCoordinates, A);
    if(temp.cartesian().dot(v.cartesian()) < 0) temp = -1 * temp;
    return temp;
}

template<int dim>
LatticeVector<dim> BiCrystal<dim>::getLatticeVectorInB(const LatticeVector<dim>& v) const
{
    VectorDimI integerCoordinates;

    if(v.lattice == &(this->B))
        return v;
    else if(v.lattice == &(this->csl))
        // V*N*csl2cslp*v
        integerCoordinates = this->matrixV() * N * csl2cslp * v;
    else
        throw(
            std::runtime_error("The input lattice vector should belong "
                               "to lattice B or the CSL"));
    auto temp = LatticeVector<dim>(integerCoordinates, B);
    if(temp.cartesian().dot(v.cartesian()) < 0) temp = -1 * temp;
    return temp;
}

template<int dim>
LatticeVector<dim> BiCrystal<dim>::getLatticeVectorInD(const LatticeVector<dim>& v) const
{
    VectorDimI integerCoordinates;

    MatrixDimI adjX = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixX());
    MatrixDimI adjV = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixV());

    if(v.lattice == &(this->A))
        // N*inv(U)*v
        integerCoordinates = N * adjX * v;
    else if(v.lattice == &(this->B))
        // M*inv(V)*v
        integerCoordinates = M * adjV * v;
    else if(v.lattice == &(this->csl))
        // N*M*csl2cslp*v
        integerCoordinates = N * M * csl2cslp * v;
    else if(v.lattice == &(this->dscl))
        return LatticeVector<dim>(v);
    else
        throw(std::runtime_error("The input lattice vector should belong to one of the four lattices of the bicrystal"));

    MatrixDimI adj_dscl2dsclp = MatrixDimIExt<IntScalarType, dim>::adjoint(dscl2dsclp);
    integerCoordinates = adj_dscl2dsclp * integerCoordinates;
    auto temp = LatticeVector<dim>(integerCoordinates, dscl);
    if(temp.cartesian().dot(v.cartesian()) < 0) temp = -1 * temp;

    // return LatticeVector<dim>(integerCoordinates,dscl);
    return temp;
}

template<int dim>
LatticeDirection<dim> BiCrystal<dim>::getLatticeDirectionInC(const LatticeVector<dim>& v) const
{
    VectorDimI integerCoordinates;

    MatrixDimI adjX = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixX());
    MatrixDimI adjV = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixV());
    MatrixDimI adjM = MatrixDimIExt<IntScalarType, dim>::adjoint(M);
    MatrixDimI adjN = MatrixDimIExt<IntScalarType, dim>::adjoint(N);
    MatrixDimI adjMN = MatrixDimIExt<IntScalarType, dim>::adjoint(M * N);
    MatrixDimI adj_csl2cslp = MatrixDimIExt<IntScalarType, dim>::adjoint(csl2cslp);

    if(v.lattice == &(this->A))
        // inv(M)*inv(U)*v
        integerCoordinates = adjM * adjX * v;
    else if(v.lattice == &(this->B))
        // inv(N)*inv(V)*v
        integerCoordinates = adjN * adjV * v;
    else if(v.lattice == &(this->csl))
        return LatticeDirection<dim>(v);
    else if(v.lattice == &(this->dscl))
        integerCoordinates = adjMN * dscl2dsclp * v;
    else
        throw(std::runtime_error(
            "The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

    integerCoordinates = adj_csl2cslp * integerCoordinates;
    auto temp = LatticeVector<dim>(integerCoordinates, csl);
    if(temp.cartesian().dot(v.cartesian()) < 0) temp = -1 * temp;
    return LatticeDirection<dim>(temp);
}
template<int dim>
LatticeDirection<dim> BiCrystal<dim>::getLatticeDirectionInD(const LatticeVector<dim>& v) const
{
    return LatticeDirection<dim>(getLatticeVectorInD(v));
}

template<int dim>
ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInA(const ReciprocalLatticeVector<dim>& rv) const
{
    VectorDimI integerCoordinates;

    MatrixDimI adjX = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixX());
    MatrixDimI adjM = MatrixDimIExt<IntScalarType, dim>::adjoint(M);
    MatrixDimI adj_csl2cslp = MatrixDimIExt<IntScalarType, dim>::adjoint(csl2cslp);
    MatrixDimI adj_dscl2dsclp = MatrixDimIExt<IntScalarType, dim>::adjoint(dscl2dsclp);

    if(rv.lattice == &(this->A))
        return ReciprocalLatticeDirection<dim>(rv);
    else if(rv.lattice == &(this->B))
        // U^-T*inverse(M)*N*V^T
        integerCoordinates = adjX.transpose() * adjM * N * (this->matrixV()).transpose() * rv;
    else if(rv.lattice == &(this->csl))
        // U^-T*inverse(M) * csl2cslp^{-T}
        integerCoordinates = adjX.transpose() * adjM * adj_csl2cslp.transpose() * rv;
    else if(rv.lattice == &(this->dscl))
        // U^-T*N*rv * dscl2dsclp^{-T}
        integerCoordinates = adjX.transpose() * N * adj_dscl2dsclp.transpose() * rv;
    else
        throw(std::runtime_error(
            "The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

    auto temp = ReciprocalLatticeVector<dim>(integerCoordinates, A);
    if(temp.cartesian().dot(rv.cartesian()) < 0) temp = -1 * temp;
    return ReciprocalLatticeDirection<dim>(temp);
}
template<int dim>
ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInB(const ReciprocalLatticeVector<dim>& rv) const
{
    VectorDimI integerCoordinates;

    MatrixDimI adjX = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixX());
    MatrixDimI adjV = MatrixDimIExt<IntScalarType, dim>::adjoint(this->matrixV());
    MatrixDimI adjM = MatrixDimIExt<IntScalarType, dim>::adjoint(M);
    MatrixDimI adjN = MatrixDimIExt<IntScalarType, dim>::adjoint(N);
    MatrixDimI adj_csl2cslp = MatrixDimIExt<IntScalarType, dim>::adjoint(csl2cslp);
    MatrixDimI adj_dscl2dsclp = MatrixDimIExt<IntScalarType, dim>::adjoint(dscl2dsclp);

    if(rv.lattice == &(this->A))
        // V^-T*inverse(N)*M*U^T
        integerCoordinates = adjV.transpose() * adjN * M * (this->matrixX()).transpose() * rv;
    else if(rv.lattice == &(this->B))
        return ReciprocalLatticeDirection<dim>(rv);
    else if(rv.lattice == &(this->csl))
        // V^-T*inverse(N)*rv*cslp2csl^T
        integerCoordinates = adjV.transpose() * adjN * adj_csl2cslp.transpose() * rv;
    else if(rv.lattice == &(this->dscl))
        // V^-T*M*rv*dsclp2dscl^T
        integerCoordinates = adjV.transpose() * M * adj_dscl2dsclp.transpose() * rv;
    else
        throw(std::runtime_error(
            "The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

    auto temp = ReciprocalLatticeVector<dim>(integerCoordinates, B);
    if(temp.cartesian().dot(rv.cartesian()) < 0) temp = -1 * temp;
    return ReciprocalLatticeDirection<dim>(temp);
}
template<int dim>
ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInC(const ReciprocalLatticeVector<dim>& rv) const
{
    VectorDimI integerCoordinates;
    MatrixDimI adj_dscl2dsclp = MatrixDimIExt<IntScalarType, dim>::adjoint(dscl2dsclp);

    if(rv.lattice == &(this->A))
        // M*U^T*rv
        integerCoordinates = M * this->matrixX().transpose() * rv;
    else if(rv.lattice == &(this->B))
        // N*V^T*rv
        integerCoordinates = N * this->matrixV().transpose() * rv;
    else if(rv.lattice == &(this->csl))
        return ReciprocalLatticeDirection<dim>(rv);
    else if(rv.lattice == &(this->dscl))
        // M*N*rv
        integerCoordinates = M * N * adj_dscl2dsclp.transpose() * rv;
    else
        throw(std::runtime_error(
            "The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

    integerCoordinates = csl2cslp.transpose() * integerCoordinates;
    auto temp = ReciprocalLatticeVector<dim>(integerCoordinates, csl);
    if(temp.cartesian().dot(rv.cartesian()) < 0) temp = -1 * temp;
    return ReciprocalLatticeDirection<dim>(temp);
}
template<int dim>
ReciprocalLatticeDirection<dim> BiCrystal<dim>::getReciprocalLatticeDirectionInD(const ReciprocalLatticeVector<dim>& rv) const
{
    VectorDimI integerCoordinates;
    MatrixDimI adj_csl2cslp = MatrixDimIExt<IntScalarType, dim>::adjoint(csl2cslp);

    MatrixDimI adjM = MatrixDimIExt<IntScalarType, dim>::adjoint(M);
    MatrixDimI adjN = MatrixDimIExt<IntScalarType, dim>::adjoint(N);

    if(rv.lattice == &(this->A))
        integerCoordinates = adjN * (this->matrixX().transpose()) * rv;
    else if(rv.lattice == &(this->B))
        integerCoordinates = adjM * (this->matrixV().transpose()) * rv;
    else if(rv.lattice == &(this->csl))
        integerCoordinates = adjN * adjM * adj_csl2cslp.transpose() * rv;
    else
        throw(std::runtime_error(
            "The input reciprocal lattice vector should belong to one of the four reciprocal lattices of the bicrystal"));

    integerCoordinates = dscl2dsclp.transpose() * integerCoordinates;
    auto temp = ReciprocalLatticeVector<dim>(integerCoordinates, dscl);
    if(temp.cartesian().dot(rv.cartesian()) < 0) temp = -1 * temp;
    return ReciprocalLatticeDirection<dim>(temp);
}

template<int dim>
LatticeVector<dim> BiCrystal<dim>::shiftTensorA(const LatticeVector<dim>& d) const
{
    if(d.lattice != &this->dscl) throw(std::runtime_error("Input vector is not a DSCL vectors"));
    return LatticeVector<dim>((LambdaA * d).eval(), *d.lattice);
}

template<int dim>
LatticeVector<dim> BiCrystal<dim>::shiftTensorB(const LatticeVector<dim>& d) const
{
    if(d.lattice != &this->dscl) throw(std::runtime_error("Input vector is not a DSCL vectors"));
    return LatticeVector<dim>((LambdaB * d).eval(), *d.lattice);
}

template<int dim>
template<typename Callback>
void BiCrystal<dim>::generateGrainBoundaries(const LatticeDirection<dim>& d, int div, GBCharacter character, Callback&& callback) const
    requires(dim == 2 || dim == 3)
{
    if(d.lattice != &A && d.lattice != &B) throw std::runtime_error("The axis does not belong to lattices A and B  ");

    constexpr IntScalarType keyScale = 1e6;
    std::optional<GBKey<keyScale>> gbKey;

    // Tilt dedup: unchanged from the original implementation. The angle-bucket test is exact here
    // because every Tilt candidate lies on a 1-parameter circle in the plane orthogonal to the axis.
    std::bitset<GBKey<keyScale>::numKeys()> seenGBs;
    seenGBs.reset();

    // Mixed dedup/key: the angle to a single reference vector is not injective over the full sphere of
    // directions (many directions share one bucket - a whole cone), so an exact integer-direction test is
    // used for deduplication, and the returned map key is a composite (angle bucket, tie-breaker) so that
    // std::map::emplace in the non-callback overload below can never silently drop two distinct Mixed GBs
    // that happen to land in the same angle bucket.
    std::set<std::array<IntScalarType, dim>> seenNormalsExact;
    constexpr IntScalarType tieBreakerMultiplier = 1'000'000;
    IntScalarType tieBreaker = 0;

    auto canonicalCoordinates = [](const Gb<dim>& gb) {
        std::array<IntScalarType, dim> coords{};
        const auto& v = gb.nA.reciprocalLatticeVector();
        for(int k = 0; k < dim; ++k) coords[k] = v(k);
        return coords;
    };

    const auto basis = d.lattice->directionOrthogonalReciprocalLatticeBasis(d, true);

    // Reciprocal direction parallel to the axis, expressed in d's own lattice (A or B, whichever the
    // caller supplied) - consistent with basis[1]/basis[2], which are also d.lattice-based. Used to build
    // Twist/Mixed candidate normals (basis[0] itself is NOT parallel to d - its only defining property is
    // basis[0].dot(d) == 1, a Bezout/duality normalization, not parallelism).
    const auto axisReciprocal = d.lattice->reciprocalLatticeDirection(d.cartesian());

    // The same axis, but pinned to lattice A specifically (regardless of whether d.lattice is A or B),
    // since Gb::nA is always expressed in A's dual. Used only by the Mixed-branch Tilt/Twist filter below,
    // which compares against gb.nA and would otherwise silently mix vectors from different lattices.
    const auto dInA = A.latticeDirection(d.cartesian());
    const auto axisReciprocalInA = A.reciprocalLatticeDirection(d.cartesian());

    if constexpr(dim == 3) {
        switch(character) {
            case GBCharacter::Tilt: {
                for(int i = -div; i <= div; ++i) {
                    for(int j = -div; j <= div; ++j) {
                        if(i == 0 && j == 0) continue;
                        const ReciprocalLatticeVector<dim> rv =
                            i * basis[1].reciprocalLatticeVector() + j * basis[2].reciprocalLatticeVector();
                        try {
                            Gb<dim> gb(*this, rv);
                            if(!gbKey) {
                                gbKey.emplace(gb);
                            }
                            const IntScalarType key = (*gbKey)(gb);
                            // Filter duplicates
                            if(seenGBs.test(key)) {
                                continue;
                            }
                            seenGBs.set(key);
                            std::invoke(callback, key, std::move(gb));
                        }
                        catch(std::runtime_error& e) {
                            Logger::warn() << e.what();
                            Logger::warn() << "Unable to form GB with normal = " << rv;
                            Logger::warn() << "moving on to next inclination";
                        }
                    }
                }
                break;
            }
            case GBCharacter::Twist: {
                // A twist boundary has no inclination freedom for a fixed misorientation: the normal is
                // forced to be parallel to the axis, i.e. axisReciprocal. div is unused for this character.
                const ReciprocalLatticeVector<dim> rv = axisReciprocal.reciprocalLatticeVector();
                try {
                    Gb<dim> gb(*this, rv);
                    constexpr IntScalarType key = 0;
                    std::invoke(callback, key, std::move(gb));
                }
                catch(std::runtime_error& e) {
                    Logger::warn() << e.what();
                    Logger::warn() << "Unable to form twist GB with normal = " << rv;
                }
                break;
            }
            case GBCharacter::Mixed: {
                for(int i = -div; i <= div; ++i) {
                    for(int j = -div; j <= div; ++j) {
                        for(int k = -div; k <= div; ++k) {
                            if(i == 0 && j == 0 && k == 0) continue;
                            const ReciprocalLatticeVector<dim> rv = i * axisReciprocal.reciprocalLatticeVector() +
                                                                     j * basis[1].reciprocalLatticeVector() +
                                                                     k * basis[2].reciprocalLatticeVector();
                            try {
                                Gb<dim> gb(*this, rv);

                                // Skip candidates that are actually pure Tilt (axis lies exactly in the GB
                                // plane: dInA.dot(nA) == 0) or pure Twist (nA parallel to the axis:
                                // nA x axisReciprocalInA == 0) - both exact integer tests. Those are already
                                // covered by the Tilt/Twist characters, so Mixed only reports genuinely
                                // mixed-character boundaries.
                                if(dInA.dot(gb.nA) == 0) continue;
                                if(gb.nA.reciprocalLatticeVector()
                                       .cross(axisReciprocalInA.reciprocalLatticeVector())
                                       .latticeVector()
                                       .squaredNorm() == 0)
                                    continue;

                                if(!gbKey) {
                                    gbKey.emplace(gb);
                                }
                                // Filter duplicates: exact test on the canonicalized (gcd-reduced) normal,
                                // since candidate normals are no longer confined to a plane.
                                if(!seenNormalsExact.insert(canonicalCoordinates(gb)).second) {
                                    continue;
                                }
                                assert(tieBreaker < tieBreakerMultiplier &&
                                       "div too large for Mixed: tie-breaker overflowed its reserved range");
                                const IntScalarType key = (*gbKey)(gb) * tieBreakerMultiplier + tieBreaker++;
                                std::invoke(callback, key, std::move(gb));
                            }
                            catch(std::runtime_error& e) {
                                Logger::warn() << e.what();
                                Logger::warn() << "Unable to form GB with normal = " << rv;
                                Logger::warn() << "moving on to next candidate";
                            }
                        }
                    }
                }
                break;
            }
        }
    }
    else if constexpr(dim == 2) {
        auto rv = basis[0].reciprocalLatticeVector();
        Gb<dim> gb = Gb<dim>(*this, rv);
        constexpr IntScalarType key = 0;
        std::invoke(callback, key, std::move(gb));
    }
}

template<int dim>
std::map<typename BiCrystal<dim>::IntScalarType, Gb<dim>> BiCrystal<dim>::generateGrainBoundaries(const LatticeDirection<dim>& d,
                                                                                                  int div,
                                                                                                  GBCharacter character) const
    requires(dim == 2 || dim == 3)
{
    std::map<IntScalarType, Gb<dim>> gbSet;
    generateGrainBoundaries(d, div, character, [&](IntScalarType key, Gb<dim>&& gb) { gbSet.emplace(key, std::move(gb)); });
    return gbSet;
}

template<int dim>
void BiCrystal<dim>::updateBoxVectors(std::array<LatticeVector<dim>, dim>& boxVectors, const double& orthogonality) const
    requires(dim == 2 || dim == 3)
{
    assert(orthogonality <= 1.0 && "The \"orthogonality\" parameter should be between 0.0 and 1.0");
    for(const auto& boxVector : boxVectors) {
        assert(&csl == boxVector.lattice && "Box vectors do not belong to the CSL.");
    }

    // Orthogonalize boxVectors[1] against the fixed boxVectors[2] (dim==3 only - in 2D there is a
    // single in-plane vector, nothing to orthogonalize against). Rather than a single one-shot integer
    // projection (which cannot always escape a tie, and only ever considers the primitive vectors
    // themselves), reuse the same "scale + RLLL search" algorithm already implemented below for
    // boxVectors[0], one dimension down: build a trivial BiCrystal<2> of the GB-plane lattice with
    // itself (Sigma=1, so its csl is exactly that lattice) purely to call its own updateBoxVectors.
    // boxVectors[2] is never modified; boxVectors[1] may grow in length as a result.
    if constexpr(dim == 3) {
        Eigen::Vector3d e1 = boxVectors[2].cartesian().normalized();
        Eigen::Vector3d e2 = (boxVectors[1].cartesian() - boxVectors[1].cartesian().dot(e1) * e1).normalized();

        // dim-1 (rather than the literal 2) is used throughout this block so that Lattice<dim-1>,
        // BiCrystal<dim-1>, etc. remain dependent on this function's own template parameter - a
        // non-dependent Lattice<2> would be looked up (and its instantiation attempted) at the point
        // BiCrystalImplementation.h is parsed, which is before Lattice.h's full definition is visible
        // (LatticeModule.h includes BiCrystal.h before Lattice.h); dim-1 defers that to the actual
        // instantiation site, by which point every header is available.
        Eigen::Matrix<double, dim - 1, dim - 1> inPlaneBasis;
        inPlaneBasis.col(0) << boxVectors[1].cartesian().dot(e1), boxVectors[1].cartesian().dot(e2);
        inPlaneBasis.col(1) << boxVectors[2].cartesian().dot(e1), boxVectors[2].cartesian().dot(e2);
        Lattice<dim - 1> inPlaneLattice(inPlaneBasis);
        BiCrystal<dim - 1> inPlaneBc(inPlaneLattice, inPlaneLattice, false);

        LatticeVector<dim - 1> v1(inPlaneBc.csl);
        v1 << 1, 0;  // == boxVectors[1], by construction
        LatticeVector<dim - 1> v2(inPlaneBc.csl);
        v2 << 0, 1;  // == boxVectors[2], fixed reference
        std::array<LatticeVector<dim - 1>, dim - 1> boxVectors2D{v1, v2};
        inPlaneBc.updateBoxVectors(boxVectors2D, orthogonality);

        // boxVectors2D[0]'s integer coordinates (a,b) are exactly the combination a*boxVectors[1] +
        // b*boxVectors[2] in the original 3D CSL - no Cartesian back-conversion needed.
        boxVectors[1] = boxVectors2D[0](0) * boxVectors[1] + boxVectors2D[0](1) * boxVectors[2];
    }

    // Adjust boxVector[0] such that it is as orthogonal as possible to boxVector[1]
    auto boxVectorTemp = boxVectors[0];
    ReciprocalLatticeDirection<dim> nC(csl);
    if constexpr(dim == 2) nC = boxVectors[1].cross();
    if constexpr(dim == 3) nC = boxVectors[1].cross(boxVectors[2]);
    auto basis = csl.planeParallelLatticeBasis(nC, true);

    int planesToExplore;
    if(abs(orthogonality) < FLT_EPSILON)
        planesToExplore = 1;
    else
        planesToExplore = nC.stacking();
    MatrixDimI boxLatticeIndices;
    boxLatticeIndices.col(0) = boxVectors[0];
    for(int i = 1; i < dim; ++i) boxLatticeIndices.col(i) = boxVectors[i] / IntegerMath<long long int>::gcd(boxVectors[i]);
    double minDotProduct = std::numbers::pi / 2;
    int minStep;

    auto boxVectorUpdated(boxVectors[0]);

    for(int i = 0; i < planesToExplore; ++i) {
        int sign = boxVectors[0].cartesian().dot(basis[0].cartesian()) > 0 ? 1 : -1;
        boxVectorTemp = boxVectors[0] + i * sign * basis[0].latticeVector();
        boxLatticeIndices.col(0) = boxVectorTemp;
        Lattice<dim> boxLattice(csl.latticeBasis * boxLatticeIndices.template cast<double>());
        ReciprocalLatticeVector<dim> rC(boxLattice);
        rC(0) = 1;
        VectorDimI temp =
            boxLatticeIndices * boxLattice.planeParallelLatticeBasis(ReciprocalLatticeDirection<dim>(rC), true)[0].latticeVector();
        boxVectorTemp = LatticeVector<dim>(temp, csl);
        double dotProduct = abs(acos(boxVectorTemp.cartesian().normalized().dot(rC.cartesian().normalized()) - FLT_EPSILON));
        if(dotProduct < minDotProduct) {
            minDotProduct = dotProduct;
            boxVectorUpdated = boxVectorTemp;
            if(dotProduct < (1 - orthogonality) * std::numbers::pi / 2) break;
        }
    }
    boxVectors[0] = boxVectorUpdated;
}

template<int dim>
std::vector<LatticeVector<dim>> BiCrystal<dim>::box(const std::array<LatticeVector<dim>, dim>& boxVectors,
                                                    const int& dsclFactor,
                                                    std::string filename,
                                                    bool orient) const
    requires(dim == 2 || dim == 3)
{
    assert(dsclFactor >= 0 && "The \"dsclFactor\" should be non-negative integer.");
    for(const auto& boxVector : boxVectors) {
        assert(&csl == boxVector.lattice && "Box vectors do not belong to the CSL.");
    }

    // Form the box lattice
    MatrixDimD C;
    for(int i = 0; i < dim; ++i) {
        C.col(i) = boxVectors[i].cartesian();
    }
    assert(abs(C.determinant()) > FLT_EPSILON && "Box volume is equal to zero.");

    ReciprocalLatticeDirection<dim> nC(csl);
    if constexpr(dim == 2) nC = boxVectors[1].cross();
    if constexpr(dim == 3) nC = boxVectors[1].cross(boxVectors[2]);

    // form the rotation matrix used to orient the system
    MatrixDimD rotation = Eigen::Matrix<double, dim, dim>::Identity();
    ;
    Eigen::Matrix<double, dim, dim - 1> orthogonalVectors;
    if(orient) {
        if constexpr(dim == 3) {
            orthogonalVectors.col(0) = C.col(1).normalized();
            orthogonalVectors.col(1) = C.col(2).normalized();
        }
        else if constexpr(dim == 2)
            orthogonalVectors.col(0) = C.col(1).normalized();

        rotation = Rotation<dim>(orthogonalVectors);
    }
    // assert((rotation*rotation.transpose()).template isApprox(Eigen::Matrix<double,dim,dim>::Identity())
    // assert((rotation*rotation.transpose()).isApprox(Eigen::Matrix<double,dim,dim>::Identity())
    assert((rotation * rotation.transpose()).isApprox(Eigen::Matrix<double, dim, dim>::Identity(), 1e-6) &&
           "Cannot orient the grain boundary. Box vectors are not orthogonal.");

    std::vector<LatticeVector<dim>> configurationA, configurationB, configurationC, configurationD;
    std::vector<LatticeVector<dim>> configuration;

    std::array<LatticeVector<dim>, dim> boxVectorsInA, boxVectorsInB, boxVectorsInD;
    // calculate boxVectors in A, B, and D
    for(int i = 0; i < dim; ++i) {
        boxVectorsInA[i] = getLatticeVectorInA(boxVectors[i]);
        boxVectorsInB[i] = getLatticeVectorInB(boxVectors[i]);
        boxVectorsInD[i] = getLatticeVectorInD(boxVectors[i]);
    }

    // prepare boxVectors for D
    auto dsclVector = getLatticeDirectionInD(boxVectors[0]).latticeVector();
    auto nD = getReciprocalLatticeDirectionInD(nC.reciprocalLatticeVector());
    if(dsclFactor != 0 && abs((dsclFactor * dsclVector).dot(nD)) < abs(boxVectorsInD[0].dot(nD)))
        boxVectorsInD[0] = dsclFactor * dsclVector;

    auto boxVectorsForA = boxVectorsInA;
    auto boxVectorsForB = boxVectorsInB;
    auto boxVectorsForC = boxVectors;
    auto boxVectorsForD = boxVectorsInD;
    boxVectorsForA[0] = 2 * boxVectorsInA[0];
    boxVectorsForB[0] = 2 * boxVectorsInB[0];
    boxVectorsForC[0] = 2 * boxVectors[0];
    boxVectorsForD[0] = 2 * boxVectorsInD[0];

    configurationA = A.box(boxVectorsForA);
    configurationB = B.box(boxVectorsForB);
    configurationC = csl.box(boxVectorsForC);
    if(dsclFactor != 0) configurationD = dscl.box(boxVectorsForD);

    LatticeVector<dim> origin(-1 * boxVectors[0]);
    for(auto& vector : configurationA) vector = vector + LatticeVector<dim>(-1 * boxVectorsInA[0]);
    for(auto& vector : configurationB) vector = vector + LatticeVector<dim>(-1 * boxVectorsInB[0]);
    for(auto& vector : configurationC) vector = vector + origin;
    for(auto& vector : configurationD) vector = vector + LatticeVector<dim>(-1 * boxVectorsInD[0]);

    configuration = configurationA;
    configuration.insert(configuration.end(), configurationB.begin(), configurationB.end());
    configuration.insert(configuration.end(), configurationC.begin(), configurationC.end());
    configuration.insert(configuration.end(), configurationD.begin(), configurationD.end());

    if(!filename.empty()) {
        std::ofstream file;
        file.open(filename);
        if(!file) std::cerr << "Unable to open file";
        file << configuration.size() << std::endl;
        file << "Lattice=\"";

        if constexpr(dim == 2) {
            file << (rotation * boxVectorsForC[0].cartesian()).transpose() << " 0 ";
            file << (rotation * boxVectorsForC[1].cartesian()).transpose() << " 0 ";
            file << " 0 0 1 ";
            file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"T T F\" origin=\"";
            file << (rotation * origin.cartesian()).transpose() << " 0.0\"" << std::endl;
            for(const auto& vector : configurationA)
                file << 1 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
            for(const auto& vector : configurationB)
                file << 2 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.05 << std::endl;
            for(const auto& vector : configurationC)
                file << 3 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.2 << std::endl;
            for(const auto& vector : configurationD)
                file << 4 << " " << (rotation * vector.cartesian()).transpose() << " " << 0.0 << "  " << 0.01 << std::endl;
        }
        else if constexpr(dim == 3) {
            file << (rotation * boxVectorsForC[0].cartesian()).transpose() << " ";
            file << (rotation * boxVectorsForC[1].cartesian()).transpose() << " ";
            file << (rotation * boxVectorsForC[2].cartesian()).transpose() << " ";
            file << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\"";
            file << (rotation * origin.cartesian()).transpose() << "\"" << std::endl;

            for(const auto& vector : configurationA)
                file << 1 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
            for(const auto& vector : configurationB)
                file << 2 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.05 << std::endl;
            for(const auto& vector : configurationC)
                file << 3 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.2 << std::endl;
            for(const auto& vector : configurationD)
                file << 4 << " " << (rotation * vector.cartesian()).transpose() << "  " << 0.01 << std::endl;
        }

        file.close();
    }
    return configuration;
}

}  // namespace oILAB

#endif
