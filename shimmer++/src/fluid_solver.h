/* This code is part of the SHIMMER project
*
* Politecnico di Torino, Dipartimento di Matematica (DISMA)
* 
* Karol Cascavita (C) 2023
* karol.cascavita@polito.it  
*/

#pragma once

#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/conservation_matrices.h"
#include "../src/assemble.h"

namespace shimmer{

void
linearized_fluid_solver(const double & tolerance, const double& dt, const double& Tm,
                        const incidence& inc, const infrastructure_graph& graph,
                        const vector_t& temp_c2_pipes, const vector_t& temp_c2_nodes,
                        vector_t& sol_time); 

} //end namespace shimmer