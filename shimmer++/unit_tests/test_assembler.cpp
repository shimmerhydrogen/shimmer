/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/conservation_matrices.h"
#include "../src/assemble.h"
#include "verify_test.h"

using triple_t = std::array<double, 3>;

using namespace shimmer;


void
make_init_graph(infrastructure_graph& igraph)
{

    std::vector<vertex_descriptor> vds;

    vds.push_back( boost::add_vertex( { "station 0", 0, 0., 0.,  0.}, igraph) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 0.,  0.}, igraph) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 0.,  0.}, igraph) );

    edge_properties ep0  = {edge_type::pipe, 0,    80000, 0.6, 9.037840034191122e-03};//0.012};
    edge_properties ep1  = {edge_type::pipe, 1,    90000, 0.6, 9.167633023370068e-03};//0.012};
    edge_properties ep2  = {edge_type::pipe, 2,   100000, 0.6, 1.113364421774635e-02};//0.012};

    /*                                            
    //           0                        *0  *1  *2              
    //         / |                     --------------         
    //        /  |                     0|  1   1                        
    //     *1/   |*0                   1| -1      -1          
    //      /____|                     2|     -1   1        
    //    1   *2   2                                         
    */

    boost::add_edge( vds[0], vds[1], ep0, igraph);
    boost::add_edge( vds[0], vds[2], ep1, igraph);
    boost::add_edge( vds[2], vds[1], ep2, igraph);
}



int main()
{

    std::vector<triple_t> ref_lhs =
                          {{{0 , 0 , 0.0009656491051934202}, 
                            {3 , 0 , 10078534.55024885},
                            {4 , 0 , 10091402.87660982}, 
                            {1 , 1 , 0.001020154052634392}, 
                            {3 , 1 ,-10078534.55024885}, 
                            {5 , 1 ,-9967287.426858675}, 
                            {2 , 2 , 0.001077080804942649}, 
                            {4 , 2 ,-10091402.87660982},                             
                            {5 , 2 , 9967287.426858675}, 
                            {0 , 3 , 1}, 
                            {1 , 3 ,-1}, 
                            {3 , 3 ,-118020405657.1786}, 
                            {0 , 4 , 1},
                            {2 , 4 ,-1},
                            {4 , 4 ,-121243672815.3844}, 
                            {1 , 5 ,-1}, 
                            {2 , 5 , 1}, 
                            {5 , 5 ,-60205728337.66899}}};

    size_t num_edges = 3;
    size_t num_vertices = 3;

    double dt = 180;

    vector_t flux (num_edges), c2_edges(num_edges);
    vector_t pressure (num_vertices), c2_vertices (num_vertices); 

    flux << 2.448496272217528e+01,
            2.171503727782473e+01,
            6.315037277824726e+00; 
    
    pressure << 5101325.0, 4977209.550248852, 4990077.876609823;


    c2_vertices << 138267.2930151191,
                138578.7460692530,
                138546.3842273756;

    c2_edges << 138422.1905008311,
                138406.1731482413,
                138562.5561693802;


    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    vector_t phi_vec = phi_vector(dt, c2_vertices, graph);
    vector_t res_friction = resistance_friction(c2_edges, flux, graph);
    vector_t res_inertia  = resistance_inertia(dt, pressure, inc, graph);
    sparse_matrix_t sAPA  = apa_matrix(c2_edges, pressure, graph, inc);

    vector_t res_vec =  res_friction + res_inertia; 
    sparse_matrix_t LHS;

    assemble_lhs(phi_vec, res_vec, sAPA, inc.matrix(), graph, LHS);
    
    std::cout << "LHS: " <<  std::endl ;
    size_t count = 0;
    for (int k = 0; k < LHS.outerSize(); ++k)
    {
        for (itor_t it(LHS,k); it; ++it, count++)
        { 
            std::cout << std::setprecision(16) << "{" << it.row() 
                      << " , " << it.col() << " , " << it.value() 
                      << "} " << std::endl ;
        }
    }

    bool lhs_pass = verify_test("LHS matrix", LHS, ref_lhs);

    return !lhs_pass; 
}