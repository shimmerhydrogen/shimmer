/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <iterator>

#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/conservation_matrices.h"
#include "../src/matrix_manipulations.h"

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;


/*  
    num_nodes        | PHI    A   Ic| | p | 
    num_pipes        | APA   -R   0 | | G |  
    num_nodes_bnd    |  Ib    0   Ib| | Gb| 
*/ 


    
sparse_matrix_t
assemble_lhs(   const vector_t & phi_vec, 
                const vector_t & res_vec,
                const sparse_matrix_t & sAPA,
                const sparse_matrix_t & sA,
                const sparse_matrix_t & sIc,
                const infrastructure_graph & graph);


vector_t
assemble_rhs(   const vector_t & rhs_continuity,
                const vector_t & rhs_momentum,
                const infrastructure_graph & graph);


void 
assemble_boundary();


std::pair<sparse_matrix_t, vector_t>
assemble(const std::pair<std::vector<triplet_t>, vector_t>& lhs_rhs_mass, 
         const std::pair<std::vector<triplet_t>, vector_t>& lhs_rhs_mom, 
         const infrastructure_graph & graph);




} //end namespace shimmer