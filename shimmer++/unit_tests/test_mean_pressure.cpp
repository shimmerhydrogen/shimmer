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
#include "verify_test.h"


using namespace shimmer;


static void
make_init_graph(infrastructure_graph& igraph)
{

 std::vector<vertex_descriptor> vds;

    vds.push_back( boost::add_vertex( { "station 0", 0, 0., 0., 100 }, igraph) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 0.,  30 }, igraph) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 0.,  60 }, igraph) );
    vds.push_back( boost::add_vertex( { "station 3", 3, 0., 0.,  80 }, igraph) );

    edge_properties ep0  = {edge_type::pipe, 0,   5, 0.7, 0.012};
    edge_properties ep1  = {edge_type::pipe, 1,   9, 0.2, 0.012};
    edge_properties ep2  = {edge_type::pipe, 2,   7, 0.3, 0.012};

    /*           0                                *0  *1  *2    
    //           |                              -------------   
    //           |*0                           0|  1           
    //           1                             1| -1  -1   1   
    //         / |                             2|         -1  
    //        /  |                             3|      1        
    //     *2/   |*1                               
    //      /    |                             
    //     2     3                                         
    */

    boost::add_edge( vds[0], vds[1], ep0, igraph);
    boost::add_edge( vds[3], vds[1], ep1, igraph);
    boost::add_edge( vds[1], vds[2], ep2, igraph);
}


int main()
{
    std::vector<double> ref_pm = {2533.333333333333,
                                    5266.666666666666,
                                    4083.333333333333};

    infrastructure_graph graph;
    make_init_graph(graph);

    vector_t pressure (num_vertices(graph)); 
    
    pressure << 2000, 3000, 5000, 7000; 

    incidence inc(graph);
    vector_t pm = average(pressure, inc);

    bool pass = verify_test(" mean pressure in pipes ", pm, ref_pm);

    return  !pass;

}