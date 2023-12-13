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

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;

class assembler
{

    /*  
        num_nodes        | PHI   Ap  I | | p | 
        num_pipes        | ADP   R   0 | | G |  
        num_nodes_bnd    |  Ib   0   Ib| | Gb| 
    */ 

    std::vector<triplet_t> t_sA; 
    std::vector<triplet_t> t_sAin;
    std::vector<triplet_t> t_sAout;
    std::vector<triplet_t> t_sPHI;
    std::vector<triplet_t> t_sADP;
    std::vector<triplet_t> t_sR;  

    incidence inc; 

    sparse_matrix_t LHS;

public:
    assembler(const infrastructure_graph& graph)
    {
        inc = incidence(graph);
        size_t num_nodes = num_vertices(graph); 
        size_t num_pipes = num_edges(graph);

        t_sA    = inc.triplets();
        t_sAin  = inc.triplets_in();
        t_sAout = inc.triplets_out();
    }


    void boundary() 
    {     
        // Bnd conditions
        // apply_offset()
        // apply_offset(t_sI, 0, num_nodes + num_pipes);
        // sparse_matrix_t sI(num_nodes, num_bnd_nodes);
        // sparse_matrix_t sIb(num_nodes, num_bnd_nodes);
    }


    void conservation_matrices(double dt, double c2,
                const vector_t& flux,
                const vector_t& pm,
                const infrastructure_graph& graph)
    {
        sparse_matrix_t sADP, sR, sPHI;
 
        t_sPHI = phi_triplets(dt, c2, graph);
        t_sADP = adp_triplets(c2, graph, inc);
        t_sR   = resistance_triplets(dt, c2, flux, pm, graph);
    }


    void assemble(const infrastructure_graph& graph)
    {
        size_t num_nodes = num_vertices(graph); 
        size_t num_pipes = num_edges(graph);

        auto apply_offset = []( std::vector<triplet_t>& triplets,
                                const std::vector<triplet_t>& vec,
                                size_t row_offset, size_t col_offset){
                                    std::transform(vec.cbegin(),vec.cend(),
                                                            std::back_inserter(triplets),
                                    [&](const triplet_t& t)
                                    { 
                                        return triplet_t(t.row() + row_offset,
                                                             t.col() + col_offset,
                                                             t.value()); 
                                    }
                                    );
                                return;
                            };

        std::vector<triplet_t> triplets; 
        triplets = t_sPHI; 
        apply_offset(triplets, t_sA,   0, num_nodes);
        apply_offset(triplets, t_sADP, num_nodes, 0);
        apply_offset(triplets, t_sR,   num_nodes, num_nodes);
        

    //-------------------------------------------------------------------------

        std::cout<< " R triplets: " <<std::endl;
        for (const auto& t : t_sR)
        {
            std::cout << std::setprecision(16) << "(" << t.row() 
                      << " , " << t.col() << " , " << t.value() 
                      << " ) " << std::endl ;
        }
        std::cout<< std::endl;
    //-------------------------------------------------------------------------

        size_t system_size = num_nodes + num_pipes; // num_bnd 
        LHS.resize(system_size, system_size);
        LHS.setFromTriplets(triplets.begin(), triplets.end());
    }


    const sparse_matrix_t& 
    incidence_in()
    { 
        return inc.matrix_in();     
    }


    const sparse_matrix_t& 
    incidence_out()
    { 
        return inc.matrix_out();     
    }


    const sparse_matrix_t& 
    LHS_matrix()
    { 
        return LHS;     
    }

};

} //end namespace shimmer