//
// Created by Nikhil Chandra Admal on 5/31/24.
//
#ifndef OILAB_PERIODICFUNCTIONIMPLEMENTATION_H
#define OILAB_PERIODICFUNCTIONIMPLEMENTATION_H
#include "FFT.h"
#include "PeriodicFunction.h"

namespace oILAB {
template <typename Scalar, int dim, int ambientDim>
PeriodicFunction<Scalar, dim, ambientDim>::PeriodicFunction(
    const Eigen::array<Eigen::Index, dim> &n,
    const Eigen::Matrix<double, ambientDim, dim> &_unitCell)
    : values(n), unitCell(_unitCell) {
  values.setZero();
    }

    template<typename Scalar, int dim, int ambientDim>
    template<typename T>
    PeriodicFunction<Scalar,dim,ambientDim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                     const Eigen::Matrix<double,ambientDim,dim>&  _unitCell,
                     const Function<T,Scalar>& fun) requires (dim == 1) :
            values(n), unitCell(_unitCell)
    {
        for (int i = 0; i < n[0]; i++)
        {
            Eigen::Vector<double,ambientDim> x=i*unitCell.col(0)/n[0];
            values(i)= fun(x);
        }
    }

    template<typename Scalar, int dim, int ambientDim>
    template<typename T>
    PeriodicFunction<Scalar,dim,ambientDim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                     const Eigen::Matrix<double,ambientDim,dim>&  _unitCell,
                     const Function<T,Scalar>& fun) requires (dim == 2) :
            values(n), unitCell(_unitCell)
    {
        for (int i = 0; i < n[0]; i++)
        {
            for (int j = 0; j < n[1]; j++)
            {
                Eigen::Vector<double,ambientDim> x=i*unitCell.col(0)/n[0] + j*unitCell.col(1)/n[1];
                values(i,j)= fun(x);
            }
        }
    }

    template<typename Scalar, int dim, int ambientDim>
    template<typename T>
    PeriodicFunction<Scalar,dim,ambientDim>::PeriodicFunction(const Eigen::array<Eigen::Index,dim>& n,
                     const Eigen::Matrix<double,ambientDim,dim>&  _unitCell,
                     const Function<T,Scalar>& fun) requires (dim == 3) :
            values(n), unitCell(_unitCell)
    {
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                for (int k = 0; k < n[2]; k++) {
                    Eigen::Vector<double,ambientDim> x= i*unitCell.col(0)/n[0] +
                                                 j*unitCell.col(1)/n[1] +
                                                 k*unitCell.col(2)/n[2];
                    values(i, j, k) = fun(x);
                }
            }
        }
    }



    template<typename Scalar, int dim, int ambientDim>
    LatticeFunction<typename PeriodicFunction<Scalar,dim,ambientDim>::dcomplex,dim,ambientDim> PeriodicFunction<Scalar,dim,ambientDim>::fft() const
    {
        Eigen::Matrix<double,ambientDim,dim> basisVectors(unitCell.transpose().completeOrthogonalDecomposition().pseudoInverse());
        LatticeFunction<dcomplex,dim,ambientDim> pfhat(values.dimensions(),basisVectors);
        FFT::fft(values.template cast<dcomplex>(),pfhat.values);

        Eigen::Matrix<double,dim,dim> unitCellGramMatrix;
        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                unitCellGramMatrix(i,j)= unitCell.col(i).dot(unitCell.col(j));
        Eigen::array<Eigen::Index,dim> n= this->values.dimensions();
        int prod= std::accumulate(std::begin(n),
                                  std::begin(n) + dim,
                                  1,
                                  std::multiplies<>{});
        pfhat.values*= (sqrt(unitCellGramMatrix.determinant())/prod);
        return pfhat;
    }

    template<typename Scalar, int dim, int ambientDim>
    double PeriodicFunction<Scalar,dim,ambientDim>::dot(const PeriodicFunction<Scalar,dim,ambientDim>& other) const
    {
        Eigen::Tensor<double,0> sum((this->values * other.values).sum());
        Eigen::Matrix<double,dim,dim> unitCellGramMatrix;
        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                unitCellGramMatrix(i,j)= unitCell.col(i).dot(unitCell.col(j));
        Eigen::array<Eigen::Index,dim> n= this->values.dimensions();
        int prod= std::accumulate(std::begin(n),
                        std::begin(n) + dim,
                        1,
                        std::multiplies<>{});
        return sum(0)*sqrt(unitCellGramMatrix.determinant())/prod;
    }

    // this needs to be corrected as we altered the definition of fft
    template<typename Scalar, int dim, int ambientDim> template <typename T>
    PeriodicFunction<Scalar,dim,ambientDim> PeriodicFunction<Scalar,dim,ambientDim>::kernelConvolution(const Function<T,Scalar>& kernel)
    {
        // pfhat - fft of the periodic function
        const auto pfhat(fft());
        Eigen::array<Eigen::Index,dim> n= values.dimensions();
        PeriodicFunction<Scalar,dim,ambientDim> output(n, unitCell);

        // fourier transform of the kernel function
        LatticeFunction<Scalar, dim, ambientDim> pkfhat(kernel.fft(n,pfhat.basisVectors));

        PeriodicFunction<dcomplex,dim,ambientDim> tempOutput(n,unitCell);
        FFT::ifft(pfhat.values*pkfhat.values,tempOutput.values);
        output.values = tempOutput.values.real();
        return output;
    }

    } // namespace oILAB
#endif // OILAB_PERIODICFUNCTIONIMPLEMENTATION_H
