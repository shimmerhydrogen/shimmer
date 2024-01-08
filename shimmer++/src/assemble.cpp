/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#include "../src/assemble.h"

namespace shimmer{


/*  
    num_nodes        | PHI    A   Ic| | p | 
    num_pipes        | APA   -R   0 | | G |  
    num_nodes_bnd    |  Ib    0   Ib| | Gb| 
*/ 


std::pair<sparse_matrix_t, vector_t>
assemble(const std::pair<std::vector<triplet_t>, vector_t>& lhs_rhs_mass, 
         const std::pair<std::vector<triplet_t>, vector_t>& lhs_rhs_mom, 
         const infrastructure_graph & graph)
{
    std::vector<triplet_t> triplets = lhs_rhs_mass.first;
    auto lhs_mom_begin = lhs_rhs_mom.first.cbegin();
    auto lhs_mom_end   = lhs_rhs_mom.first.cend();
    triplets.insert(triplets.begin(),lhs_mom_begin, lhs_mom_end);

    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);
    size_t system_size = num_nodes + num_pipes;// + num_pipes_ext; // num_bnd 

    sparse_matrix_t LHS(system_size, system_size);    
    LHS.setFromTriplets(triplets.begin(), triplets.end()); 
    
    vector_t rhs = vector_t::Zero(system_size);
    rhs.head(num_nodes) =  lhs_rhs_mass.second;   
    rhs.segment(num_nodes, num_pipes) = lhs_rhs_mom.second;   

    return std::make_pair(LHS, rhs);
}


void assemble_boundary() 
{     
    // Bnd conditions
    // apply_offset()
    // apply_offset(t_sI, 0, num_nodes + num_pipes);
    // sparse_matrix_t sI(num_nodes, num_bnd_nodes);
    // sparse_matrix_t sIb(num_nodes, num_bnd_nodes);
}




} //end namespace shimmer