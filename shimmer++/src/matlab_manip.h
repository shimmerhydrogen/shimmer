 /* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it  
 */

#pragma once

#include <Eigen/Dense>
#include "../src/incidence_matrix.h"
#include "infrastructure_graph.h"

namespace shimmer{

using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

matrix_t
build_x_nodes(const infrastructure_graph& g);


} // namespace shimmer