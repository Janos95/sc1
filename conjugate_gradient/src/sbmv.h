#pragma once

#include <Eigen/Sparse>
#include <Eigen/Core>


extern "C" 
{
#include <mkl.h>
}

namespace sc1
{

template<typename DerivedA, typename DerivedB>
DerivedB mklSbmv(const DerivedA& A, const  DerivedB& x)
{
    DerivedB y(x.rows());
    
    if constexpr(std::is_same<typename DerivedA::Scalar, float>::value)
        cblas_ssbmv (CblasRowMajor, CblasUpper, A.rows(), 1, 1.0f, A.valuePtr(), 3, x.data(), 1, 0.0f, y.data(), 1);
      
    else if constexpr(std::is_same<typename DerivedA::Scalar, double>::value)
        cblas_dsbmv (CblasRowMajor, CblasUpper, A.rows(), 1, 1.0, A.valuePtr(), 3, x.data(), 1, 0.0, y.data(), 1);
    
    return y;
}

template<typename DerivedA, typename DerivedB>
DerivedB eigenSbmv(const DerivedA& A, const  DerivedB& x)
{
    return A * x;
}

}
