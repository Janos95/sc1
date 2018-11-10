#pragma once

//#include <boost/timer/timer.hpp>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>


namespace sc1
{
    
template<typename DerivedA>
class ConjugateGradient
{
public:        
        template<typename DerivedB, typename MatVecMult>
        DerivedB solve(const Eigen::MatrixBase<DerivedB>& b, const MatVecMult mult)
        {
            //boost::timer::auto_cpu_timer t("Running cg took %w seconds\n");
            using namespace Eigen;
            
            MatrixBase<DerivedB> x = DerivedB::Zero(b.rows());
            MatrixBase<DerivedB> r = b - mult(m_A, x);
            MatrixBase<DerivedB> p = r;
            typename DerivedA::Scalar rsold = r.dot(r);
            
            for(int i = 0; i < m_iterations; ++i)
            {
                DerivedB Ap = mult(m_A, p);
                typename DerivedA::Scalar alpha = rsold / p.dot(Ap);
                x +=  alpha * p;
                r -= alpha * Ap;
                typename DerivedA::Scalar rsnew = r.dot(r);
                
                if(std::sqrt(rsnew) < m_minResidual)
                    break;
                
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
        
private:
	Eigen::MatrixBase<DerivedA> m_A;
	int m_iterations = 100;
        typename DerivedA::Scalar m_minResidual = typename DerivedA::Scalar(1e-10);
};

}


//class Base<T>

//class Derived : Base<Derived>
