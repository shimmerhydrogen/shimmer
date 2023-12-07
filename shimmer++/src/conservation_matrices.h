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
#include "infrastructure_graph.h"
#include "../src/geometry_properties.h"
#include "../src/incidence_matrix.h"

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


void
average(const  vector_t& pressure, const sparse_matrix_t& incidence_in,
        const sparse_matrix_t& incidence_out,  vector_t& pm);

void
phi_matrix(const double & dt, const double& c2, const infrastructure_graph& g, sparse_matrix_t& mat);


void 
adp_matrix(const double & c2, const infrastructure_graph& g,
                const sparse_matrix_t& incidence_in,
                const sparse_matrix_t& incidence_out,
                sparse_matrix_t& mat);


void
resistance_matrix(const double & dt, const double& c2,
                  const vector_t & flux,
                  const vector_t & mean_pressure,
                  const infrastructure_graph  & g,
                  sparse_matrix_t& mat );

} //end namespace shimmer
