/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#pragma once

#include <Eigen/Sparse>
#include "infrastructure_graph.h"

template<typename T> 
using sparse_matrix_t = Eigen::SparseMatrix<T>; 


template<typename T>
sparse_matrix_t<T>
incidence_matrix(const infrastructure_graph& g)
{
    using triplet_t = Eigen::Triplet<T>;

    std::vector<triplet_t> triplets;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
        auto pipe = g[*itor];   
        auto u = source(*itor, g);
        auto v = target(*itor, g);
        triplets.push_back(triplet_t(g[u].node_num, pipe.branch_num, T(1)));
        triplets.push_back(triplet_t(g[v].node_num,   pipe.branch_num, T(-1)));
    }

    sparse_matrix_t<T> mat(num_vertices(g), num_edges(g));
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    return mat;
}
