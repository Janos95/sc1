#include "conjugate_gradient.h"
#include "build_problem.h"

#include <Eigen/Sparse>

#include <boost/timer/timer.hpp>

#include <cstdlib>


using Scalar = double;

int main(int argc, char** argv)
{
    boost::timer::auto_cpu_timer t("Whole pipeline took %w seconds\n");
    
    if(argc!=2)
    {
        std::cerr << "Error: expected one and only one argument.\n";
        return -1;
    }	
    
    int n = std::atoi(argv[1]);
    auto [A, b] = sc1::buildPoissonProblem(n);

    sc1::ConjugateGradient cg;
    cg.setMatrix(A);
    
    auto x = cg.solve(b, mklSbmv);
    
    return 0;
}
