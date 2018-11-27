#pragma once

#include "build_problem.h"

#include <boost/timer/timer.hpp>

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace sc1
{
    
template<typename MatrixT>
class ConjugateGradient
{
public:        
        template<typename VectorT, typename MatVecMult>
        VectorT solve(const VectorT& b, const MatVecMult mult)
        {
            boost::timer::auto_cpu_timer t("Running cg took %w seconds\n");
            
            using namespace Eigen;
            
            if(!m_iterations)
                m_iterations = b.rows();
            
            m_iterationsUsed = m_iterations;
            
            VectorT x = VectorT::Zero(b.rows());
            VectorT r = b - mult(m_A, x);
            VectorT p = r;
            typename MatrixT::Scalar rsold = r.dot(r);
            
            for(int i = 0; i < m_iterations; ++i)
            {
                VectorT Ap = mult(m_A, p);
                typename MatrixT::Scalar alpha = rsold / p.dot(Ap);
                x +=  alpha * p;
                r -= alpha * Ap;
                typename MatrixT::Scalar rsnew = r.dot(r);
                
                if(std::sqrt(rsnew) < m_minResidual)
                {
                    m_iterationsUsed = i+1;
                    fmt::printf("Terminating cg after %d iterations because residual %.16g is smaller than %.16g\n", i+1, std::sqrt(rsnew), m_minResidual);
                    break;
                }
                
                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
            }
            
            return x;
        }
        
        void setMatrix(const MatrixT& A)
        {
            m_A = A;
        }
        
        void setIterations(int iterations)
        {
            m_iterations = iterations;
        }
        
        void setMinResidual(typename MatrixT::Scalar minResidual)
        {
            m_minResidual = minResidual;
        }
        
        int iterations()
        {
            return m_iterationsUsed;
        }
        
private:
	MatrixT m_A;
	int m_iterations = 0;
        int m_iterationsUsed;
        typename MatrixT::Scalar m_minResidual = typename MatrixT::Scalar(1e-10);
};

}

