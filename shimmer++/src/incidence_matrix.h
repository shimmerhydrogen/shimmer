/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#pragma once

#include <Eigen/Sparse>
#include "infrastructure_graph.h"

using triplet_t = Eigen::Triplet<double>;
using sparse_matrix_t = Eigen::SparseMatrix<double>; 

sparse_matrix_t
incidence_matrix(const infrastructure_graph& g){


    std::vector<triplet_t> triplets;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
        auto pipe = g[*itor]; 
        std::cout << pipe << std::endl;
        triplets.push_back(triplet_t(pipe.from, pipe.branch_num, 1.0));
        triplets.push_back(triplet_t(pipe.to,   pipe.branch_num, -1.0));
    }

    std::cout << num_vertices(g) << std::endl;
    std::cout << num_edges(g)    << std::endl;

    sparse_matrix_t mat(num_vertices(g), num_edges(g));
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    return mat;
}
