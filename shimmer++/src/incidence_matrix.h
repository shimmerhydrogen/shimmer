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
    incidence_in_triplets(const infrastructure_graph& g);

    void
    incidence_out_triplets(const infrastructure_graph& g);

    void
    incidence_triplets(const infrastructure_graph& g);

    void 
    incidence_matrix_in(const infrastructure_graph& g);

    void
    incidence_matrix_out(const infrastructure_graph& g);

    void
    incidence_matrix(const infrastructure_graph& g);

public:
    incidence(const infrastructure_graph& g)
    {
        incidence_triplets(g);
        incidence_matrix(g);
    }

    incidence(){}

    size_t triplets_size();    
    size_t triplets_in_size(); 
    size_t triplets_out_size();
    size_t triplets_size()  const; 
    size_t triplets_in_size() const; 
    size_t triplets_out_size()const; 

    const sparse_matrix_t& matrix();     
    const sparse_matrix_t& matrix_in();   
    const sparse_matrix_t& matrix_out();   
    const sparse_matrix_t& matrix() const;      
    const sparse_matrix_t& matrix_in()  const;   
    const sparse_matrix_t& matrix_out() const;  

    const std::vector<triplet_t>&  triplets() const;     
    const std::vector<triplet_t>&  triplets_in() const; 
    const std::vector<triplet_t>&  triplets_out() const; 

    const std::vector<triplet_t>& triplets();
    const std::vector<triplet_t>& triplets_in();
    const std::vector<triplet_t>& triplets_out();  
    
};



} //end namespace shimmer