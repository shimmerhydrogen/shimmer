 /* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#pragma once

#include <iostream>
#include <Eigen/Sparse>

namespace shimmer{


using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using itor_t = sparse_matrix_t::InnerIterator;


std::vector<triplet_t>
build_triplets(const vector_t& vec, size_t row_offset = 0,
                                    size_t col_offset = 0);


sparse_matrix_t
build_matrix(const vector_t& vec);


std::vector<triplet_t>
build_triplets(const sparse_matrix_t& sMat, size_t row_offset = 0,
                                            size_t col_offset = 0);

}// end namespace shimmer
