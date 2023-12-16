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
    num_nodes        | PHI    A   I | | p | 
    num_pipes        | APA   -R   0 | | G |  
    num_nodes_bnd    |  Ib    0   Ib| | Gb| 
*/ 


    
void assemble_lhs(  const vector_t& phi_vec, 
                    const vector_t& res_vec,
                    const sparse_matrix_t & sAPA,
                    const sparse_matrix_t & sA,
                    const infrastructure_graph& graph,
                    sparse_matrix_t& LHS)
{
    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);

    auto t_sPHI = build_triplets( phi_vec);
    auto t_sR   = build_triplets(-res_vec, num_nodes, num_nodes);
    auto t_sAPA = build_triplets( sAPA,    num_nodes, 0);
    auto t_sA   = build_triplets( sA,      0, num_nodes );

    std::vector<triplet_t> triplets =  t_sPHI; 
    triplets.insert(triplets.begin(), t_sAPA.begin(), t_sAPA.end());
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());

    size_t system_size = num_nodes + num_pipes; // num_bnd 
    LHS.resize(system_size, system_size);
    LHS.setFromTriplets(triplets.begin(), triplets.end());
}

/*
void assemble_rhs(flux, pressure)
{

    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);

    phi_matrix(dt, c2, graph, PHI);
    resistance_matrix(dt, c2, flux, pm, graph, R);
    resistance_inertia_matrix(dt, c2, flux, pm, graph, Ri);

    rhs.head(num_nodes) = PHI * pressure;
	rhs.segment(num_nodes, num_pipes)  = R * flux + Ri * (flux_n - flux); // -(Rf.*abs(G_k(:,k)).*G_k(:,k) + Ri.*G_n);
}
*/
void assemble_boundary() 
{     
    // Bnd conditions
    // apply_offset()
    // apply_offset(t_sI, 0, num_nodes + num_pipes);
    // sparse_matrix_t sI(num_nodes, num_bnd_nodes);
    // sparse_matrix_t sIb(num_nodes, num_bnd_nodes);
}




} //end namespace shimmer