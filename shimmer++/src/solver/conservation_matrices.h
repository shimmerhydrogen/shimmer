/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <Eigen/Sparse>
#include "infrastructure_graph.h"
#include "solver/geometry_properties.h"
#include "solver/incidence_matrix.h"
#include "solver/pipe_calculator.h"
#include "solver/matlab_manip.h"
#include "solver/gas_law.h"
#include "solver/variable.h"

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
