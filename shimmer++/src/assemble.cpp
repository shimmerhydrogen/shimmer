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


    
sparse_matrix_t
assemble_lhs(   const vector_t & phi_vec, 
                const vector_t & res_vec,
                const sparse_matrix_t & sAPA,
                const sparse_matrix_t & sA,
                const sparse_matrix_t & sIc,
                const infrastructure_graph & graph)
{
    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);
    size_t num_pipes_ext = num_pipes;

    auto t_sPHI = build_triplets( phi_vec);
    auto t_sR   = build_triplets(-res_vec, num_nodes, num_nodes);
    auto t_sAPA = build_triplets( sAPA,    num_nodes, 0);
    auto t_sA   = build_triplets( sA,      0, num_nodes );
    auto t_sIc  = build_triplets( sIc,     0, num_nodes + num_pipes);

    std::vector<triplet_t> triplets =  t_sPHI; 
    triplets.insert(triplets.begin(), t_sAPA.begin(), t_sAPA.end());
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());
    //triplets.insert(triplets.begin(), t_sIc.begin(), t_sIc.end());

    size_t system_size = num_nodes + num_pipes + num_pipes_ext; // num_bnd 
    sparse_matrix_t LHS(system_size, system_size);
    LHS.setFromTriplets(triplets.begin(), triplets.end());
    return LHS;
}


vector_t
assemble_rhs(   const vector_t & rhs_continuity,
                const vector_t & rhs_momentum,
                const infrastructure_graph & graph)
{
    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);
    size_t system_size = num_nodes + num_pipes;// + num_pipes_ext; // num_bnd 

    vector_t rhs = vector_t::Zero(system_size);
    rhs.head(num_nodes) =  rhs_continuity;   
    rhs.segment(num_nodes, num_pipes) = rhs_momentum;   
    return rhs;
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