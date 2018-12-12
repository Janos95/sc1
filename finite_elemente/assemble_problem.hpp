//
// Created by janos und karsten GmbH Co. KG on 12/1/18. This is free software but please dont kill anybody!
// Use it in responsible way!
//

#pragma once

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <range/v3/all.hpp>

#include <fmt/core.h>

#include <cmath>



namespace sc1
{

template<typename T>
auto assembleProblem(const int N)
{
    //boost::timer::auto_cpu_timer t("Building the poisson problem took %w seconds\n");
    using namespace ranges;

    auto mapToIndex = [=](const int i, const int j){ return j * (N + 1) + i;};

    Eigen::Matrix<T, Eigen::Dynamic, 1> b((2*N+1)*(N+1), 1);


    for(const auto& [m,n]: view::cartesian_product(view::ints(0,2*N+1), view::ints(0,N+1)))
    {

        bool onBoundary = m == 0 || m == 2*N || n == 0 || n == N;
        b(mapToIndex(m,n), 0) = onBoundary ? T{1} : T{2};
    }

    //set corners correctly
    b(0,0) = T{1}/T{2};
    b(2*N*(N+1),0) = T{1}/T{2};
    b(N,0) = T{1}/T{2};
    b(2*N*(N+1)+N,0) = T{1}/T{2};


    return b;
}


};