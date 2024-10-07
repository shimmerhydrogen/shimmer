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
#include "../src/pipe_calculator.h"
#include "../src/matlab_manip.h"
#include "../src/gas_law.h"
#include "MATLAB_GERG_functions.hpp"
#include "../src/variable.h"

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;



vector_t
average(const  vector_t& pressure, const incidence& inc);


vector_t
phi_vector( const double & dt, const vector_t& c2,
            const infrastructure_graph& g);


sparse_matrix_t
phi_matrix( const double & dt, const vector_t& c2,
            const infrastructure_graph& g);


sparse_matrix_t
adp_matrix(const vector_t& c2, const infrastructure_graph& g,
            const incidence& inc);


vector_t
resistance_inertia( const double & dt, const vector_t & pressure,
                    const incidence& inc, const infrastructure_graph & g);


vector_t
resistance_friction(const double& temperature,
                    const vector_t& mu,
                    const vector_t& c2,
                    const vector_t & flux,
                    const infrastructure_graph & g);




} //end namespace shimmer
