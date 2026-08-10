//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_LATTICEFUNCTION_H
#define OILAB_LATTICEFUNCTION_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

// Lattice function defined on the lattice points
namespace oILAB {

template <typename Derived, typename Scalar> class Function;

// NOTE: the default for ambientDim lives on the class definition in
// PeriodicFunction.h -- repeating it here would redefine it whenever both
// headers are included.
template <typename Scalar, int dim, int ambientDim> class PeriodicFunction;

template <typename Scalar, int dim, int ambientDim = dim> class LatticeFunction {
  using dcomplex = std::complex<double>;

public:
  const Eigen::Matrix<double, ambientDim, dim> basisVectors;
  Eigen::Tensor<Scalar, dim> values;
  explicit LatticeFunction(
      const Eigen::array<Eigen::Index, dim> &n,
      const Eigen::Matrix<double, ambientDim, dim> &_basisVectors);

  template <typename T>
  LatticeFunction(
      const Eigen::array<Eigen::Index, dim> &n,
      const Eigen::Matrix<double, ambientDim, dim> &_basisVectors,
      const Function<T, Scalar> &fun) requires (dim == 1);

  template <typename T>
  LatticeFunction(
      const Eigen::array<Eigen::Index, dim> &n,
      const Eigen::Matrix<double, ambientDim, dim> &_basisVectors,
      const Function<T, Scalar> &fun) requires (dim == 2);

  template <typename T>
  LatticeFunction(
      const Eigen::array<Eigen::Index, dim> &n,
      const Eigen::Matrix<double, ambientDim, dim> &_basisVectors,
      const Function<T, Scalar> &fun) requires (dim == 3);

  std::complex<double>
  dot(const LatticeFunction<std::complex<double>, dim, ambientDim> &other) const;

  PeriodicFunction<dcomplex, dim, ambientDim> ifft() const;
    };

    template<typename Scalar, int dim, int ambientDim>
    LatticeFunction<Scalar, dim, ambientDim> operator*(const LatticeFunction<Scalar,dim,ambientDim>& lf1, const LatticeFunction<Scalar,dim,ambientDim>& lf2);
    } // namespace oILAB

#include "LatticeFunctionImplementation.h"

#endif //OILAB_LATTICEFUNCTION_H
