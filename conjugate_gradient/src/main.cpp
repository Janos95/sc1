#define _SILENCE_CXX17_NEGATORS_DEPRECATION_WARNING

#include "conjugate_gradient.h"
#include "build_problem.h"
#include "sbmv.h"
#include "mkl_tridiagonal_solver.h"

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <pcl/visualization/pcl_plotter.h>
#include <vtk-6.2/vtkChart.h>

#include <boost/timer/timer.hpp>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <functional>


using Scalar = double; //not be integral

int global;

Scalar computeCgError(Scalar nx, int iterations)
{
    int n = nx;
    
    auto [A, b] = sc1::buildPoissonProblem<Scalar>(n);
    sc1::ConjugateGradient<decltype(A)> cg;
    cg.setIterations(global);
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
        std::cerr << "Error: expected one and only one argument.\n";
        return -1;
    }	
    
    int n = std::atoi(argv[1]); 
    std::istringstream os(argv[2]); 
    Scalar minResidual;
    os >> minResidual;
 

    pcl::visualization::PCLPlotter * plotter = new pcl::visualization::PCLPlotter ();

    std::vector<char> red = {(char)255,(char)0,(char)0,(char)255};
    std::vector<char> blue = {(char)0,(char)0,(char)255,(char)255};
    
    plotter->addPlotData (computeTriDiagError, 100, 500, "TriDiagSolver", 100, vtkChart::POINTS, blue);
    
    for(const int iterations: {1, 5, 500})
    {
        std::ostringstream ss;
        ss << "CgSolver with " << iterations << " iterations";
        
        using namespace std::placeholders;
        std::function<double(double)> cgErrorIterations = std::bind(computeCgError, _1, iterations);
        
        plotter->addPlotData (cgErrorIterations, 100, 500, ss.str().c_str(), 100, vtkChart::POINTS, red);
    }
    
    plotter->plot ();
    
    delete plotter;
    
    return 0;
}
