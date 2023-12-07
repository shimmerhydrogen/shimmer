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

using sparse_matrix_t = Eigen::SparseMatrix<double>; 


sparse_matrix_t
incidence_matrix_out(const infrastructure_graph& g);


sparse_matrix_t
incidence_matrix_in(const infrastructure_graph& g);


sparse_matrix_t
incidence_matrix(const infrastructure_graph& g);