#include "conjugate_gradient.h"
#include "build_problem.h"
#include "sbmv.h"

// #include <boost/timer/timer.hpp>

#include <cstdlib>
#include <iostream>


using Scalar = double; //not be integral

int main(int argc, char** argv)
{
    //boost::timer::auto_cpu_timer t("Whole pipeline took %w seconds\n");
    
    if(argc!=2)
    {
        std::cerr << "Error: expected one and only one argument.\n";
        return -1;
    }	
    
    int n = std::atoi(argv[1]);
    
    
    auto [A, b] = sc1::buildPoissonProblem<Scalar>(n); //TODO: why cannot be const?
    const auto solution = sc1::computeExactSolution<Scalar>(n);
    
   // std::cout << A << std::endl;
    //std::cout << b << std::endl;

    sc1::ConjugateGradient<decltype(A)> cg;
    cg.setMatrix(A);
    
    auto x = cg.solve(b, sc1::eigenSbmv<decltype(A), decltype(b)>);
    
//     std::cout << "computed solution" << x.head<10>() << std::endl;
//    std::cout << "exact solution" << b << std::endl;
    
    std::cout << "The error in the maximum norm is " << (x-solution).lpNorm<Eigen::Infinity>() << std::endl;
    return 0;
}
