#pragma once

#include <boost/timer/timer.hpp>

extern "C" 
{
#include <mkl.h>
}

namespace sc1
{

template<typename DerivedA, typename DerivedB>
DerivedB mklTriDiagonalSolver(const DerivedA& A, const  DerivedB& b)
{
    boost::timer::auto_cpu_timer t("MKL tridiagonal solver took %w seconds\n");
    
    DerivedB diag = A.diagonal();
    DerivedB upperDiag = A.block(0, 1, A.rows() - 1, A.cols() - 1).eval().diagonal();
    DerivedB lowerDiag = A.block(1, 0, A.rows() - 1, A.cols() - 1).eval().diagonal();
    
//     std::cout << diag << std::endl << std::endl;
//     std::cout << upperDiag << std::endl << std::endl;
//     std::cout << lowerDiag <<  std::endl << std::endl;
    DerivedB x = b;
 
    if constexpr(std::is_same<typename DerivedA::Scalar, float>::value)
        auto i = LAPACKE_sgtsv(CblasColMajor , A.rows() , 1 ,  lowerDiag.data(), diag.data() ,  upperDiag.data(), x.data() , b.rows());
      
    else if constexpr(std::is_same<typename DerivedA::Scalar, double>::value)
        auto i = LAPACKE_dgtsv(CblasColMajor , A.rows() , 1 ,  lowerDiag.data(), diag.data() ,  upperDiag.data(), x.data() , b.rows());
    
    return x;
}

}
