//
// Created by Nikhil Chandra Admal on 7/4/23.
//

#ifndef OILAB_PERIODICFUNCTION_H
#define OILAB_PERIODICFUNCTION_H

#include "Eigen/Dense"
#include "Function.h"
#include "unsupported/Eigen/CXX11/Tensor"
#include <iomanip>

namespace oILAB {

template <typename Scalar, int dim, int ambientDim = dim> class LatticeFunction;

template <typename Derived, typename Scalar> class Function;

template <typename Scalar, int dim, int ambientDim = dim> class PeriodicFunction {
public:
  using dcomplex = std::complex<double>;
  const Eigen::Matrix<double, ambientDim, dim> unitCell;
  Eigen::Tensor<Scalar, dim> values;

  // Construct a zero periodic function
  explicit PeriodicFunction(
      const Eigen::array<Eigen::Index, dim> &n,
      const Eigen::Matrix<double, ambientDim, dim> &_unitCell);

  // generates a periodic function from a function centered at the center of a
  // lattice
  template <typename T>
  PeriodicFunction(const Eigen::array<Eigen::Index, dim> &n,
                   const Eigen::Matrix<double, ambientDim, dim> &_unitCell,
                   const Function<T, Scalar> &fun) requires (dim == 1);

  template <typename T>
  PeriodicFunction(const Eigen::array<Eigen::Index, dim> &n,
                   const Eigen::Matrix<double, ambientDim, dim> &_unitCell,
                   const Function<T, Scalar> &fun) requires (dim == 2);

  template <typename T>
  PeriodicFunction(const Eigen::array<Eigen::Index, dim> &n,
                   const Eigen::Matrix<double, ambientDim, dim> &_unitCell,
                   const Function<T, Scalar> &fun) requires (dim == 3);

  LatticeFunction<dcomplex, dim, ambientDim> fft() const;

  double dot(const PeriodicFunction<Scalar, dim, ambientDim> &other) const;

  template <typename T>
  PeriodicFunction<Scalar, dim, ambientDim>
  kernelConvolution(const Function<T, Scalar> &kernel);
    };

    template<typename Scalar, int dim, int ambientDim>
    std::basic_ostream<char>& operator<<(std::basic_ostream<char>& s, const PeriodicFunction<Scalar, dim, ambientDim>& fun) requires (dim==2)
    {
        auto n = fun.values.dimensions();
        assert(n.size() == dim);
        s << n[0]*n[1] << std::endl;
        s << std::endl;
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                Eigen::Vector<double, ambientDim> x = i * fun.unitCell.col(0) / n[0] +
                                                      j * fun.unitCell.col(1) / n[1];
                const Eigen::IOFormat fmt(15, 0, " ", "", " ", "", "", "");
                s << x.transpose().format(fmt) << std::setw(25) << std::setprecision(15) << fun.values(i, j)
                  << std::endl;
            }
        }
        return s;
    }

    template<typename Scalar, int dim, int ambientDim, typename T>
    std::basic_ostream<char>& operator<<(std::basic_ostream<char>& s, const PeriodicFunction<Scalar, dim, ambientDim>& fun) requires (dim==3)
    {
        auto n = fun.values.dimensions();
        assert(n.size() == dim);
        s << n[0]*n[1]*n[2] << std::endl;
        s << std::endl;
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                for (int k = 0; k < n[2]; k++) {
                    Eigen::Vector<double, ambientDim> x = i * fun.unitCell.col(0) / n[0] +
                                                          j * fun.unitCell.col(1) / n[1] +
                                                          k * fun.unitCell.col(2) / n[2];
                    const Eigen::IOFormat fmt(15, 0, " ", " ", " ", "", "", "");
                    s << x.transpose().format(fmt) << std::setw(25) << std::setprecision(15) << fun.values(i,j,k)
                      << std::endl;
                }
            }
        }
        return s;
    }

    } // namespace oILAB
#include "PeriodicFunctionImplementation.h"

#endif //OILAB_PERIODICFUNCTION_H
