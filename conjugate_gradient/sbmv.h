#pragma once

#include <boost/timer/timer.hpp>

namespace sc1
{

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> mklSbmv(const Eigen::SparseMatrix<T>& A, const  Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    if constexpr(std::is_same<T, float>::value
        cblas_ssbmv (CblasColMajor, 'u', n, 1, 1.0f, const float a, const MKL_INT lda, const float *x, const MKL_INT incx, const float beta, float *y, const MKL_INT incy);
      
    else if constexpr(std::is_same<T, double>::value)
        cblas_dsbmv (CblasColMajor, 'u', n, 1, 1.0, const float a, const MKL_INT lda, const float *x, const MKL_INT incx, const float beta, float *y, const MKL_INT incy);
    
    return 
}
	
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> eigenSbmv(const Eigen::SparseMatrix<T>& A, const  Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    return A * x;
}
