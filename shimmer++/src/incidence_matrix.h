/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#pragma once

#include <Eigen/Sparse>
#include "infrastructure_graph.h"

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;


class incidence
{
    sparse_matrix_t mat_; 
    sparse_matrix_t mat_in_;
    sparse_matrix_t mat_out_; 

    std::vector<triplet_t> triplets_;
    std::vector<triplet_t> triplets_in_;
    std::vector<triplet_t> triplets_out_;


    void
    compute_triplets(const infrastructure_graph& g);

    void
    compute_matrix(const infrastructure_graph& g);

public:
    incidence(const infrastructure_graph& g)
    {
        compute_triplets(g);
        compute_matrix(g);
    }

    incidence(){}

    const sparse_matrix_t& matrix();     
    const sparse_matrix_t& matrix_in();   
    const sparse_matrix_t& matrix_out();   
    const sparse_matrix_t& matrix() const;      
    const sparse_matrix_t& matrix_in()  const;   
    const sparse_matrix_t& matrix_out() const;  

};



} //end namespace shimmer