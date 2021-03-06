#pragma once

#include <cmath>

#include <boost/timer/timer.hpp>

#include <Eigen/Sparse>
#include <Eigen/Core>

namespace sc1
{

template<typename T>
auto buildPoissonProblem(const int n)
{	
    boost::timer::auto_cpu_timer t("Building the poisson problem took %w seconds\n");
    
    Eigen::SparseMatrix<T> A(n-2,n-2);
    Eigen::Matrix<T, Eigen::Dynamic, 1> b(n-2);
    {
        
        std::vector<Eigen::Triplet<T>> coefficients; 
        coefficients.reserve(3*n-2);
        
        for(int i = 0; i < n-2; ++i)
        {
            b(i) = std::sin(2 * M_PI * T(i+1) / (n-1));
            
            coefficients.emplace_back(i, i, 2);
            
            if(i > 0)
            {
                coefficients.emplace_back(i, i-1, -1); // lower
                coefficients.emplace_back(i-1, i, -1); // upper
            }
        }
        
        A.setFromTriplets(coefficients.begin(), coefficients.end());
        A *= T(1.0) * ((n-1)*(n-1));
        A.makeCompressed();
    }
    
    return std::make_pair(std::move(A), std::move(b));
}

template<typename T>
auto computeExactSolution(const int n) 
{	
    boost::timer::auto_cpu_timer t("Computing the exact solution took %w seconds\n");

    Eigen::Matrix<T, Eigen::Dynamic, 1> solution(n-2);

    for(int i = 0; i < n-2; ++i)
    {
        solution(i) = 1/(4*M_PI*M_PI) * std::sin(2 * M_PI * T(i+1) / (n-1));
    }

    return solution;
}

}
