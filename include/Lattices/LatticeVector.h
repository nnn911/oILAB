/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */

#ifndef gbLAB_LatticeVector_h_
#define gbLAB_LatticeVector_h_

#include <array>
#include "LatticeModule.h"

namespace oILAB {
/*! \brief LatticeVector class
 *
 *  The LatticeVector<dim> class describes a lattice vector in a lattice
 * */
template <int dim>
class LatticeVector
    : public Eigen::Matrix<typename LatticeCore<dim>::IntScalarType, dim, 1> {
  typedef Eigen::Matrix<typename LatticeCore<dim>::IntScalarType, dim, 1>
      BaseType;
  BaseType &base();
  const BaseType &base() const;

public:
  //        static constexpr double roundTol=FLT_EPSILON;
  typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
  typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
  typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
  typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
  typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

  const Lattice<dim> *lattice;

  LatticeVector();
  LatticeVector(const Lattice<dim> &lat);
  LatticeVector(const VectorDimD &d, const Lattice<dim> &lat);
  LatticeVector(const VectorDimI &d, const Lattice<dim> &lat);
  LatticeVector(const LatticeVector<dim> &other) = default;
  LatticeVector(LatticeVector<dim> &&other) = default;

  LatticeVector<dim> &operator=(const LatticeVector<dim> &other);
  LatticeVector<dim> &operator=(LatticeVector<dim> &&other);
  LatticeVector<dim> operator+(const LatticeVector<dim> &other) const;
  LatticeVector<dim> &operator+=(const LatticeVector<dim> &other);
  LatticeVector<dim> operator-(const LatticeVector<dim> &other) const;
  LatticeVector<dim> &operator-=(const LatticeVector<dim> &other);
  LatticeVector<dim> operator*(const IntScalarType &scalar) const;

  IntScalarType dot(const ReciprocalLatticeVector<dim> &other) const;
  IntScalarType dot(const ReciprocalLatticeDirection<dim> &other) const;
  VectorDimD cartesian() const;

  ReciprocalLatticeDirection<dim>
  cross(const LatticeVector<dim> &other) const requires (dim == 2) {
    assert(lattice == other.lattice &&
           "LatticeVectors belong to different Lattices.");
    return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(
        (VectorDimI() << 0, 0).finished(), *lattice));
  }
  ReciprocalLatticeDirection<dim>
  cross(const LatticeVector<dim> &other) const requires (dim == 3) {
    assert(lattice == other.lattice &&
           "LatticeVectors belong to different Lattices.");
    return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(
        static_cast<VectorDimI>(*this).cross(static_cast<VectorDimI>(other)),
        *lattice));
  }

  ReciprocalLatticeDirection<dim>
  cross() const requires (dim == 2) {
    return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(
        (VectorDimI() << -(*this)(1), (*this)(0)).finished(), *lattice));
  }
  ReciprocalLatticeDirection<dim>
  cross() const requires (dim == 3) {
    return ReciprocalLatticeDirection<dim>(ReciprocalLatticeVector<dim>(
        (VectorDimI() << -(*this)(1), (*this)(0), 0).finished(), *lattice));
  }

  static void modulo(
      LatticeVector<dim> &input, const std::array<LatticeVector<dim>, dim> &basis,
      const VectorDimD &shift = VectorDimD::Zero()) requires (dim == 2 || dim == 3) {
    if constexpr (dim == 3) {
      double det = (basis[0].cross(basis[1])).dot(basis[2]);
      assert(abs(det) > FLT_EPSILON);
      auto normal = basis[1].cross(basis[2]);
      input = input - floor((double)input.dot(normal) / basis[0].dot(normal) - shift(0)) * basis[0];
      assert((double)(input.dot(normal)) / basis[0].dot(normal) <= shift(0) + 1.0 &&
             (double)(input.dot(normal)) / basis[0].dot(normal) >= shift(0));
      normal = basis[2].cross(basis[0]);
      input = input - floor((double)input.dot(normal) / basis[1].dot(normal) - shift(1)) * basis[1];
      assert((double)(input.dot(normal)) / basis[1].dot(normal) <= shift(1) + 1.0 &&
             (double)(input.dot(normal)) / basis[1].dot(normal) >= shift(1));
      normal = basis[0].cross(basis[1]);
      input = input - floor((double)input.dot(normal) / basis[2].dot(normal) - shift(2)) * basis[2];
      assert((double)(input.dot(normal)) / basis[2].dot(normal) <= shift(2) + 1.0 &&
             (double)(input.dot(normal)) / basis[2].dot(normal) >= shift(2));
    } else {
      auto normal = basis[1].cross();
      input = input - input.dot(normal) / basis[0].dot(normal) * basis[0];
      normal = basis[0].cross();
      input = input - input.dot(normal) / basis[1].dot(normal) * basis[1];
    }
  }

  static void modulo(
      VectorDimD &input, const std::array<LatticeVector<dim>, dim> &basis,
      const VectorDimD &shift = VectorDimD::Zero()) requires (dim == 2 || dim == 3) {
    if constexpr (dim == 3) {
      Eigen::Matrix3d L;
      L.col(0) = basis[0].cartesian();
      L.col(1) = basis[1].cartesian();
      L.col(2) = basis[2].cartesian();
      Eigen::Vector3d inputCoordinates = ((L.inverse() * input).array() - shift.array()).floor();
      input = input - L * inputCoordinates;
    }
  }

    };

    template<int dim>
    LatticeVector<dim> operator*(const typename LatticeVector<dim>::IntScalarType& scalar, const LatticeVector<dim>& L);

    template<int dim>
    LatticeVector<dim> operator*(const int& scalar, const LatticeVector<dim>& L);
    } // namespace oILAB

#include "LatticeVectorImplementation.h"

#endif
