#pragma once

#include "build_problem.h"

#include <boost/timer/timer.hpp>

def EIGEN_USE_MKL_ALL

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace sc1
{
    
template<typename DerivedA>
class ConjugateGradient
{
public:        
        template<typename DerivedB, typename MatVecMult>
        DerivedB solve(const DerivedB& b, const MatVecMult mult)
        {
            boost::timer::auto_cpu_timer t("Running cg took %w seconds\n");
            
            using namespace Eigen;
            
            if(!m_iterations)
                m_iterations = b.rows();
            
            m_iterationsUsed = m_iterations;
            
            DerivedB x = DerivedB::Zero(b.rows());
            DerivedB r = b - mult(m_A, x);
            DerivedB p = r;
            typename DerivedA::Scalar rsold = r.dot(r);
            
            for(int i = 0; i < m_iterations; ++i)
            {
                DerivedB Ap = mult(m_A, p);
                typename DerivedA::Scalar alpha = rsold / p.dot(Ap);
                x +=  alpha * p;
                r -= alpha * Ap;
                typename DerivedA::Scalar rsnew = r.dot(r);
                
                if(std::sqrt(rsnew) < m_minResidual)
                {
                    m_iterationsUsed = i+1;
                    printf("Terminating cg after %d iterations because residual %.16g is smaller than %.16g\n", i+1, std::sqrt(rsnew), m_minResidual);
                    break;
                }
                
                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
            }
            
            return x;
        }
        
        void setMatrix(const DerivedA& A)
        {
            m_A = A;
        }
        
        void setIterations(int iterations)
        {
            m_iterations = iterations;
        }
        
        void setMinResidual(typename DerivedA::Scalar minResidual)
        {
            m_minResidual = minResidual;
        }
        
        int iterations()
        {
            return m_iterationsUsed;
        }
        
private:
	DerivedA m_A;
	int m_iterations = 0;
        int m_iterationsUsed;
        typename DerivedA::Scalar m_minResidual = typename DerivedA::Scalar(1e-10);
};

}

