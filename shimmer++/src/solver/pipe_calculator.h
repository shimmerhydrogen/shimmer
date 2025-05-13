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

#include <unordered_map>
#include <Eigen/Sparse>
#include "../src/infrastructure_graph.h"

namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


double
viscosity(const double& Temperature,
          const vector_t & X);


double
friction_factor_average(const edge_properties& pipe ,const double & Temperature,
                        const double & flux, const double & mu);


double 
inertia_resistance( const edge_properties& pipe, const double& dt,
                    const double& mean_pressure) ;

double 
friction_resistance(const  edge_properties& pipe,
                    double temperature,
                    double mu,
                    double c2, 
                    double flux);


} //end namespace shimmer
