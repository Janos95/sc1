#pragma once

#include <boost/timer/timer.hpp>

extern "C" 
{
#include <mkl.h>
}

namespace sc1
{

template<typename MatrixT, typename VectorT>
VectorT mklTriDiagonalSolver(const MatrixT& A, const  VectorT& b)
{
    boost::timer::auto_cpu_timer t("MKL tridiagonal solver took %w seconds\n");
    
    VectorT diag = A.diagonal();
    VectorT upperDiag = A.block(0, 1, A.rows() - 1, A.cols() - 1).eval().diagonal();
    VectorT lowerDiag = A.block(1, 0, A.rows() - 1, A.cols() - 1).eval().diagonal();
    
//     std::cout << diag << std::endl << std::endl;
//     std::cout << upperDiag << std::endl << std::endl;
//     std::cout << lowerDiag <<  std::endl << std::endl;
    VectorT x = b;
 
    if constexpr(std::is_same<typename MatrixT::Scalar, float>::value)
        auto i = LAPACKE_sgtsv(CblasColMajor , A.rows() , 1 ,  lowerDiag.data(), diag.data() ,  upperDiag.data(), x.data() , b.rows());
      
    else if constexpr(std::is_same<typename MatrixT::Scalar, double>::value)
        auto i = LAPACKE_dgtsv(CblasColMajor , A.rows() , 1 ,  lowerDiag.data(), diag.data() ,  upperDiag.data(), x.data() , b.rows());
    
    return x;
}

}
