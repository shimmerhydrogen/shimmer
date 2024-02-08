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
resistance_friction(const double& temperature, const vector_t& c2,
                    const vector_t & flux,
                    const infrastructure_graph & g);

using pair_trip_vec_t = std::pair<std::vector<triplet_t>, vector_t>;


pair_trip_vec_t
momentum(const double& dt, const double& temperature,
         const vector_t& flux, const vector_t& flux_old,
         const vector_t& pressure_nodes, const vector_t& pressure_pipes,
         const vector_t& c2, const incidence& inc,
         const infrastructure_graph & graph);


pair_trip_vec_t
continuity(const double& dt, const double& temperature,
        const vector_t& pressure, 
        const vector_t& pressure_old,
        const vector_t& c2,
        const incidence& inc,
        const infrastructure_graph & graph);


pair_trip_vec_t
boundary(size_t num_nodes_, size_t num_pipes,
        const vector_t& p_in,
        const vector_t& flux,
        const vector_t& flux_ext,
        const vector_t& vel,
        const vector_t& inlet_nodes);

} //end namespace shimmer
