/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_ReciprocalLatticeVector_cpp_
#define gbLAB_ReciprocalLatticeVector_cpp_

#include "../../include/IO/Logger.h"
#include "../../include/Lattices/LatticeModule.h"

namespace oILAB {

template <int dim>
typename ReciprocalLatticeVector<dim>::BaseType &
ReciprocalLatticeVector<dim>::base() {
  return *this;
    }

    template <int dim>
    const typename ReciprocalLatticeVector<dim>::BaseType& ReciprocalLatticeVector<dim>::base() const
    {
        return *this;
    }
    
    template <int dim>
    ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const Lattice<dim> &lat) :
    /* init */ BaseType(VectorDimI::Zero()),
    /* init */ lattice(lat)
    { }

    template <int dim>
    ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const VectorDimD &d, const Lattice<dim> &lat) :
    /* init */ BaseType(LatticeCore<dim>::integerCoordinates(d, lat.latticeBasis.transpose())),
    /* init */ lattice(lat)
    { }

    template <int dim>
    ReciprocalLatticeVector<dim>::ReciprocalLatticeVector(const VectorDimI& other, const Lattice<dim> &lat) :
    /* init base */ BaseType(other),
    /* init base */ lattice(lat)
    {}

    template <int dim>
    ReciprocalLatticeVector<dim>& ReciprocalLatticeVector<dim>::operator=(const ReciprocalLatticeVector<dim> &other)
    {
        assert(&lattice == &other.lattice && "ReciprocalLatticeVector<dim> belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVector<dim>& ReciprocalLatticeVector<dim>::operator=(ReciprocalLatticeVector<dim> &&other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() = other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVector<dim> ReciprocalLatticeVector<dim>::operator+(const ReciprocalLatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        VectorDimI temp= static_cast<VectorDimI>(*this) + static_cast<VectorDimI>(other);
        return ReciprocalLatticeVector<dim>(temp, lattice);
    }

    template <int dim>
    ReciprocalLatticeVector<dim>& ReciprocalLatticeVector<dim>::operator+=(const ReciprocalLatticeVector<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() += other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVector<dim> ReciprocalLatticeVector<dim>::operator-(const ReciprocalLatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        VectorDimI temp= static_cast<VectorDimI>(*this) - static_cast<VectorDimI>(other);
        return ReciprocalLatticeVector<dim>(temp, lattice);
    }

    template <int dim>
    ReciprocalLatticeVector<dim>& ReciprocalLatticeVector<dim>::operator-=(const ReciprocalLatticeVector<dim> &other)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        base() -= other.base();
        return *this;
    }

    template <int dim>
    ReciprocalLatticeVector<dim> ReciprocalLatticeVector<dim>::operator*(const IntScalarType& scalar) const
    {
        VectorDimI temp= static_cast<VectorDimI>(*this) * scalar;
        return ReciprocalLatticeVector<dim>(temp, lattice);
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::IntScalarType ReciprocalLatticeVector<dim>::dot(const LatticeVector<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return static_cast<VectorDimI>(*this).dot(static_cast<VectorDimI>(other));
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::IntScalarType ReciprocalLatticeVector<dim>::dot(const LatticeDirection<dim> &other) const
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return dot(other.latticeVector());
    }

    

    template <int dim>
    typename ReciprocalLatticeVector<dim>::VectorDimD ReciprocalLatticeVector<dim>::cartesian() const
    {
        return lattice.reciprocalBasis * this->template cast<double>();
    }
    

    template <int dim>
    typename ReciprocalLatticeVector<dim>::IntScalarType ReciprocalLatticeVector<dim>::closestPlaneIndexOfPoint(const VectorDimD &P) const
    {
      if (this->squaredNorm() == 0) {
        Logger::debug() << "A null ReciprocalLatticeVector cannot be "
                           "used to compute planeIndexOfPoint";
        assert(false && "A null ReciprocalLatticeVector cannot be "
                        "used to compute planeIndexOfPoint");
      }
      const double hd(cartesian().dot(P));
      return std::lround(hd);
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::IntScalarType ReciprocalLatticeVector<dim>::planeIndexOfPoint(const VectorDimD &P) const
    {
        assert(this->squaredNorm() > 0 && "A null ReciprocalLatticeVector cannot be used to compute planeIndexOfPoint");
        const double hd(cartesian().dot(P));
        const IntScalarType h(std::lround(hd));
        if (fabs(hd - h) > FLT_EPSILON)
        {
          Logger::debug() << "P=" << P.transpose();
          Logger::debug() << "r=" << this->cartesian().transpose();
          Logger::debug() << "hd=" << std::setprecision(15) << std::scientific
                          << hd;
          Logger::debug() << "h=" << h;
          Logger::debug() << "P in not on a lattice plane.";
          assert(false && "P in not on a lattice plane.");
        }
        return h;
    }

    template <int dim>
    typename ReciprocalLatticeVector<dim>::IntScalarType ReciprocalLatticeVector<dim>::planeIndexOfPoint(const LatticeVector<dim> &P) const
    {
      if (this->squaredNorm() == 0) {
        Logger::debug() << "A null ReciprocalLatticeVector cannot be "
                           "used to compute planeIndexOfPoint";
        assert(false && "A null ReciprocalLatticeVector cannot be "
                        "used to compute planeIndexOfPoint");
      }
      return dot(P);
    }


    //Operator starts here
    template <int dim>
    ReciprocalLatticeVector<dim> operator*(const typename ReciprocalLatticeVector<dim>::IntScalarType& scalar, const ReciprocalLatticeVector<dim> &L)
    {
        return L * scalar;
    }
    template <int dim>
    ReciprocalLatticeVector<dim> operator*(const int& scalar, const ReciprocalLatticeVector<dim> &L)
    {
        return L * scalar;
    }

    template<int dim>
    LatticeDirection<dim>
    ReciprocalLatticeVector<dim>::cross(const ReciprocalLatticeVector<dim>& other) const requires (dim==2)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeDirection<dim>(
                LatticeVector<dim>((VectorDimI() << 0,0).finished(), lattice));
    }

    template<int dim>
    LatticeDirection<dim>
    ReciprocalLatticeVector<dim>::cross() const requires (dim==2) {
        return LatticeDirection<dim>(
                LatticeVector<dim>((VectorDimI() << -(*this)(1), (*this)(0)).finished(), lattice));
    }

    template<int dim>
    LatticeDirection<dim>
    ReciprocalLatticeVector<dim>::cross(const ReciprocalLatticeVector<dim>& other) const requires (dim==3)
    {
        assert(&lattice == &other.lattice && "LatticeVectors belong to different Lattices.");
        return LatticeDirection<dim>(LatticeVector<dim>(static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)), lattice));
    }
    template<int dim>
    LatticeDirection<dim>
    ReciprocalLatticeVector<dim>::cross() const requires (dim==3) {
        return LatticeDirection<dim>(
                LatticeVector<dim>((VectorDimI() << -(*this)(1), (*this)(0), 0).finished(), lattice));
    }

template class ReciprocalLatticeVector<1>;
template ReciprocalLatticeVector<1> operator*(const typename ReciprocalLatticeVector<1>::IntScalarType& scalar, const ReciprocalLatticeVector<1> &L);
template ReciprocalLatticeVector<1> operator*(const int& scalar, const ReciprocalLatticeVector<1> &L);
template class ReciprocalLatticeVector<2>;
template ReciprocalLatticeVector<2> operator*(const typename ReciprocalLatticeVector<2>::IntScalarType&scalar, const ReciprocalLatticeVector<2> &L);
template ReciprocalLatticeVector<2> operator*(const int& scalar, const ReciprocalLatticeVector<2> &L);

template class ReciprocalLatticeVector<3>;
template ReciprocalLatticeVector<3> operator*(const typename ReciprocalLatticeVector<3>::IntScalarType& scalar, const ReciprocalLatticeVector<3> &L);
template ReciprocalLatticeVector<3> operator*(const int& scalar, const ReciprocalLatticeVector<3> &L);

template class ReciprocalLatticeVector<4>;
template ReciprocalLatticeVector<4> operator*(const typename ReciprocalLatticeVector<4>::IntScalarType& scalar, const ReciprocalLatticeVector<4> &L);
template ReciprocalLatticeVector<4> operator*(const int& scalar, const ReciprocalLatticeVector<4> &L);
template class ReciprocalLatticeVector<5>;
template ReciprocalLatticeVector<5> operator*(const typename ReciprocalLatticeVector<5>::IntScalarType& scalar, const ReciprocalLatticeVector<5> &L);
template ReciprocalLatticeVector<5> operator*(const int& scalar, const ReciprocalLatticeVector<5> &L);
} // namespace oILAB
#endif
