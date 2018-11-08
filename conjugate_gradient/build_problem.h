#pragma once

namespace sc1
{

template<typename T>
auto buildPoissonProblem(const int n)
{	
    auto_cpu_timer t("Building the poisson problem took %w seconds\n");
    
    Eigen::SparseMatrix<T> A;
    {
      std::vector<Eigen::Triplet<T>> coefficients; 
      coefficients.reserve(3*n-2);
      
      for(int i = 0; i < n; ++i)
      {
        coefficients.emplace_back(i, i, 2);
        
        if(i > 0)
        {
          coefficients.emplace_back(i, i-1, -1);
          coefficients.emplace_back(i-1, i, -1);
        }
      }
     
      A.setFromTriplets(coefficients.begin(), coefficients.end());
      A *= T(1.0) / (n-1)^2;
    }
    
    Eigen::VectorXd b;     
    
    return std::make_pair(std::move(A), std::move(b));
}

}
