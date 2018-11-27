#pragma once

#include <Eigen/Sparse>
#include <Eigen/Core>


extern "C" 
{
#include <mkl.h>
}

namespace sc1
{

template<typename MatrixT, typename VectorT>
VectorT mklSbmv(const MatrixT& A, const  VectorT& x)
{
    VectorT y(x.rows());
    
    if constexpr(std::is_same<typename MatrixT::Scalar, float>::value)
        cblas_ssbmv (CblasRowMajor, CblasUpper, A.rows(), 1, 1.0f, A.valuePtr(), 3, x.data(), 1, 0.0f, y.data(), 1);
      
    else if constexpr(std::is_same<typename MatrixT::Scalar, double>::value)
        cblas_dsbmv (CblasRowMajor, CblasUpper, A.rows(), 1, 1.0, A.valuePtr(), 3, x.data(), 1, 0.0, y.data(), 1);
    
    return y;
}

template<typename MatrixT, typename VectorT>
VectorT eigenSbmv(const MatrixT& A, const  VectorT& x)
{
    return A * x;
}

}
