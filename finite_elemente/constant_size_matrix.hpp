/*
 * Copyright (c) 2018, Karsten Evers and Janos Meny
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the
 * above copyright notice, this list of conditions
 * and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the
 * names of its contributors may be used to endorse or
 * promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
 * AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */


#pragma once

#include <range/v3/all.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class ConstantSizeMatrix;
using Eigen::SparseMatrix;


namespace Eigen {
    namespace internal {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<>
        struct traits<ConstantSizeMatrix> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
        {};
    }
}



// only square matrices possible
class ConstantSizeMatrix : public Eigen::EigenBase<ConstantSizeMatrix> {
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    template<typename Rhs>
    Eigen::Product<ConstantSizeMatrix, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const {
        return Eigen::Product<ConstantSizeMatrix, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:

    const StorageIndex n;

    const Scalar equal;
    const Scalar neighbor; 

    ConstantSizeMatrix(const StorageIndex n) :
            n(n),
            equal(4.0), // einmal
            neighbor(-1.0) // 4 mal
    {}

    StorageIndex rows() const { return (2*n-1)*(n-1); }

    StorageIndex cols() const { return (2*n-1)*(n-1); }

};


namespace Eigen {
    namespace internal {
        template<typename Rhs>
        struct generic_product_impl<ConstantSizeMatrix, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
                : generic_product_impl_base<ConstantSizeMatrix,Rhs,generic_product_impl<ConstantSizeMatrix,Rhs> >
        {
            typedef typename Product<ConstantSizeMatrix,Rhs>::Scalar Scalar;
            typedef typename Product<ConstantSizeMatrix,Rhs>::StorageIndex StorageIndex;

            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const ConstantSizeMatrix& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha==Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                const StorageIndex n = lhs.n;
                const auto mapToIndex = [=](const StorageIndex x, const StorageIndex y){ return x  + y * (n - 1);};
                const Scalar equal = lhs.equal;
                const Scalar neighbor = lhs.neighbor;


#pragma omp parallel for
                for(StorageIndex x = 1; x < n-2; ++x)
                {
                    for(StorageIndex y = 1; y < 2*n-2; ++y)
                    {
                        StorageIndex i = mapToIndex(x,y);

                        dst[i] += neighbor * rhs[i+1] //right neighbor
                                  + neighbor * rhs[i-1] // left neighbor
                                  + neighbor * rhs[i+n-1] // upper neighbor
                                  + neighbor * rhs[i-(n-1)] // lower neighbor
/*                                  + neighbor * rhs[i+n-1+1] // upper right neighbor
                                  + neighbor * rhs[i-(n-1)-1] // lower left neighbor*/
/*                                  + neighbor * rhs[i-(n-1)+1] // lower right neighbor*/
             /*                     + neighbor * rhs[i+(n-1)-1] // upper left neighbor*/
                                  + equal * rhs[i];
                    }
                }


                //edges
                {
                    //upper edge
                    for(StorageIndex i = mapToIndex(1,2*n-2); i <= mapToIndex(n-3,2*n-2); ++i)
                    {

                        dst[i] += neighbor * rhs[i+1] //right neighbor
                                  + neighbor * rhs[i-1] // left neighbor
                                  + neighbor * rhs[i-(n-1)] // lower neighbor
                          /*        + neighbor * rhs[i-(n-1)-1] // lower left neighbor*/
                       /*           + neighbor * rhs[i-(n-1)+1] // lower right neighbor*/
                                  + equal * rhs[i];
                    }

                    //lower edge
                    for(StorageIndex i = mapToIndex(1,0); i <= mapToIndex(n-3,0); ++i)
                    {
                        dst[i] += neighbor * rhs[i+1] //right neighbor
                                  + neighbor * rhs[i-1] // left neighbor
                                  + neighbor * rhs[i+n-1] // upper neighbor
                               /*   + neighbor * rhs[i+n-1+1] // upper right neighbor*/
                                /*  + neighbor * rhs[i+(n-1)-1] // upper left neighbor*/
                                  + equal * rhs[i];
                    }

                    //right edge
                    for(StorageIndex i = mapToIndex(n-2,1); i <= mapToIndex(n-2,2*n-3); i+=n-1)
                    {
                        dst[i] += neighbor * rhs[i-1] // left neighbor
                                  + neighbor * rhs[i+n-1] // upper neighbor
                                  + neighbor * rhs[i-(n-1)] // lower neighbor
                       /*           + neighbor * rhs[i-(n-1)-1] // lower left neighbor*/
                          /*        + neighbor * rhs[i+(n-1)-1] // upper left neighbor*/
                                  + equal * rhs[i];
                    }

                    //left edge
                    for(StorageIndex i = mapToIndex(0,1); i <= mapToIndex(0,2*n-3); i+=n-1)
                    {
                        dst[i] += neighbor * rhs[i+1] //right neighbor
                                  + neighbor * rhs[i+n-1] // upper neighbor
                                  + neighbor * rhs[i-(n-1)] // lower neighbor
                              /*    + neighbor * rhs[i+n-1+1] // upper right neighbor*/
                               /*   + neighbor * rhs[i-(n-1)+1] // lower right neighbor*/
                                  + equal * rhs[i];
                    }


                }


                //corners
                {
                    //upper left coner
                    StorageIndex upperLeftStorageIndex = mapToIndex(0,2*n-2);
                    dst[upperLeftStorageIndex] += equal * rhs[upperLeftStorageIndex]
                    + neighbor * rhs[upperLeftStorageIndex + 1]
                    + neighbor * rhs[upperLeftStorageIndex-(n-1)]
                   /* + neighbor * rhs[upperLeftStorageIndex-(n-1)+1]*/;

                    //upper right corner
                    StorageIndex upperRightStorageIndex = mapToIndex(n-2,2*n-2);
                    dst[upperRightStorageIndex] += equal * rhs[upperRightStorageIndex]
                    + neighbor * rhs[upperRightStorageIndex - 1]
                    + neighbor * rhs[upperRightStorageIndex-(n-1)]
                 /*   + neighbor * rhs[upperRightStorageIndex-(n-1)-1]*/;

                    //lower right corner
                    StorageIndex lowerRightStorageIndex = mapToIndex(n-2,0);
                    dst[lowerRightStorageIndex] += equal * rhs[lowerRightStorageIndex]
                    + neighbor * rhs[lowerRightStorageIndex - 1]
                    + neighbor * rhs[lowerRightStorageIndex+(n-1)]
                    /*+ neighbor * rhs[lowerRightStorageIndex+(n-1)-1]*/;

                    //lower left corner
                    StorageIndex lowerLeftStorageIndex = mapToIndex(0,0);
                    dst[lowerLeftStorageIndex] += equal * rhs[lowerLeftStorageIndex]
                    + neighbor * rhs[lowerLeftStorageIndex + 1]
                    + neighbor * rhs[lowerLeftStorageIndex+(n-1)]
                  /*  + neighbor * rhs[lowerLeftStorageIndex+(n-1)+1]*/;
                }

            }
        };
    }
}

