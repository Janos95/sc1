#pragma once


namespace sc1
{
    
template<typename MatrixType>
class ConjugateGradient
{
public:
        void setMatrix(const MatrixType& A)
        {
            m_A = A;
        }
        
        template<typename VectorType, typename Mult>
        VectorType solve(const VectorType& b, const Mult mult)
        {
            boost::timer::auto_cpu_timer t("Running cg took %w seconds\n");
            
            Eigen::VectorXd x(b.rows());
            Eigen::ConjugateGradient<MatrixType, Lower|Upper> cg;
            cg.compute(m_A);
            x = cg.solve(b);
            std::cout << "estimated error: " << cg.error()      << std::endl;
            return x;
        }
        
private:
	MatrixType m_A
	int m_numIters;
};

}
