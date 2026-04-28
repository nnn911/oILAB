/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */

#ifndef gbLAB_ReciprocalLatticeVectorBase_h_
#define gbLAB_ReciprocalLatticeVectorBase_h_

#include "LatticeModule.h"

namespace oILAB {
template <int dim>
class ReciprocalLatticeVector
    : public Eigen::Matrix<typename LatticeCore<dim>::IntScalarType, dim, 1> {
  typedef Eigen::Matrix<typename LatticeCore<dim>::IntScalarType, dim, 1>
      BaseType;
  BaseType &base();
  const BaseType &base() const;

public:
  typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
  typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
  typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
  typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
  typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

  const Lattice<dim> &lattice;

  ReciprocalLatticeVector(const Lattice<dim> &lat);
  ReciprocalLatticeVector(const VectorDimD &d, const Lattice<dim> &lat);
  ReciprocalLatticeVector(const VectorDimI &d, const Lattice<dim> &lat);
  ReciprocalLatticeVector(const ReciprocalLatticeVector<dim> &other) = default;
  ReciprocalLatticeVector(ReciprocalLatticeVector<dim> &&other) = default;

  ReciprocalLatticeVector<dim> &
  operator=(const ReciprocalLatticeVector<dim> &other);
  ReciprocalLatticeVector<dim> &operator=(ReciprocalLatticeVector<dim> &&other);
  ReciprocalLatticeVector<dim>
  operator+(const ReciprocalLatticeVector<dim> &other) const;
  ReciprocalLatticeVector<dim> &
  operator+=(const ReciprocalLatticeVector<dim> &other);
  ReciprocalLatticeVector<dim>
  operator-(const ReciprocalLatticeVector<dim> &other) const;
  ReciprocalLatticeVector<dim> &
  operator-=(const ReciprocalLatticeVector<dim> &other);
  ReciprocalLatticeVector<dim> operator*(const IntScalarType &scalar) const;
  IntScalarType dot(const LatticeVector<dim> &other) const;
  IntScalarType dot(const LatticeDirection<dim> &other) const;
  VectorDimD cartesian() const;
  IntScalarType closestPlaneIndexOfPoint(const VectorDimD &P) const;
  IntScalarType planeIndexOfPoint(const VectorDimD &P) const;
  IntScalarType planeIndexOfPoint(const LatticeVector<dim> &P) const;

  LatticeDirection<dim>
  cross(const ReciprocalLatticeVector<dim> &other) const requires (dim == 2);

  LatticeDirection<dim> cross() const requires (dim == 2);

  LatticeDirection<dim>
  cross(const ReciprocalLatticeVector<dim> &other) const requires (dim == 3);

  LatticeDirection<dim> cross() const requires (dim == 3);
    };
    
    template<int dim>
    ReciprocalLatticeVector<dim> operator*(const typename ReciprocalLatticeVector<dim>::IntScalarType& scalar, const ReciprocalLatticeVector<dim>& L);
    template<int dim>
    ReciprocalLatticeVector<dim> operator*(const int& scalar, const ReciprocalLatticeVector<dim>& L);

    } // namespace oILAB

#include "ReciprocalLatticeVectorImplementation.h"
#endif
