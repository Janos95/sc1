// #def EIGEN_USE_MKL_ALL


#include "conjugate_gradient.h"
#include "build_problem.h"
#include "sbmv.h"
#include "mkl_tridiagonal_solver.h"

#include <Eigen/Sparse>
#include <Eigen/Core>

// #include <pcl/visualization/pcl_plotter.h>
// #include <vtk-6.2/vtkChart.h>

#include <fmt/format.h>

#include <boost/timer/timer.hpp>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <functional>


using namespace ranges;
using namespace std::placeholders;

using Scalar = double; 

Scalar computeCgError(int iterations, Scalar nx)
{
    int n = nx;
    
    auto [A, b] = sc1::buildPoissonProblem<Scalar>(n);
    sc1::ConjugateGradient<decltype(A)> cg;
    cg.setIterations(iterations);
    cg.setMatrix(A);
    
    auto x = cg.solve(b, sc1::mklSbmv<decltype(A), decltype(b)>);
    auto solution = sc1::computeExactSolution<Scalar>(n);
    
    return (x-solution).lpNorm<Eigen::Infinity>();
}

Scalar computeTriDiagError(Scalar nx)
{
    int n = nx;
    
    auto [A, b] = sc1::buildPoissonProblem<Scalar>(n);
    
    auto x = sc1::mklTriDiagonalSolver(A,b);
    auto solution =sc1:: computeExactSolution<Scalar>(n);
    
    return (x-solution).lpNorm<Eigen::Infinity>();
}

int main(int argc, char** argv)
{
    boost::timer::auto_cpu_timer t("Whole pipeline took %w seconds\n");
    
    if(argc!=3)
    {
        std::cerr << "Error: expected 2 and only 2 argument.\n";
        return -1;
    }	
    
    int n = std::atoi(argv[1]); 
    std::istringstream os(argv[2]); 
    Scalar minResidual;
    os >> minResidual;
 
    Plotter plotter;
    plotter.setXAxis(view::linear_distribute(5.f, n , 100)); // problem sizes
    
    for(const auto& [iterations, color]: {{1, red}, {5, green}, {500, blue}})
    {
        std::string plotName = fmt::format("CgSolver with {} iterations", iterations);
        const auto function = std::bind(computeCgError, iterations, _1);
        plotter.plot (function, color, plotName);
    }
    
    plotter.spin();
    
    
//     auto plotter = std::make_unique<pcl::visualization::PCLPlotter>();
// 
//     std::vector<char> red = {(char)255,(char)0,(char)0,(char)255};
//     std::vector<char> blue = {(char)0,(char)0,(char)255,(char)255};
//     
//     plotter->addPlotData (computeTriDiagError, 100, 500, "TriDiagSolver", 100, vtkChart::POINTS, blue);
//     
//     for(const int iterations: {1, 5, 500})
//     {
// //         std::ostringstream ss;
// //         ss << "CgSolver with " << iterations << " iterations";
//         std::string plotName = fmt::format("CgSolver with {} iterations", iterations);
//         
//         global = iterations;
//         plotter->addPlotData (computeCgError, 100, 500, plotName.c_str(), 100, vtkChart::POINTS, red);
//     }
//     
//     plotter->plot ();
    
}
