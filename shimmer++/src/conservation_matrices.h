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


auto
phi_triplets(const double & dt, const double& c2, const infrastructure_graph& g);


void
phi_matrix(const double & dt, const double& c2, const infrastructure_graph& g,
           sparse_matrix_t& mat);


auto 
adp_triplets(const double & c2, const infrastructure_graph& g,
             const incidence& inc);


void
adp_matrix(const double & c2, const infrastructure_graph& g,
            const incidence& inc,
            sparse_matrix_t& mat);


auto
resistance_triplets(const double & dt, const double& c2,
                  const vector_t & flux,
                  const vector_t & mean_pressure,
                  const infrastructure_graph  & g);


void
resistance_matrix(const double & dt, const double& c2,
                  const vector_t & flux,
                  const vector_t & mean_pressure,
                  const infrastructure_graph  & g,
                  sparse_matrix_t& mat );

} //end namespace shimmer
