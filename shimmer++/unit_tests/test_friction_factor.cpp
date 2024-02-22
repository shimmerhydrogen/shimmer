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
#include "../src/pipe_calculator.h"
#include "verify_test.h"
#include "../src/viscosity.h"

using triple_t = std::array<double, 3>;
using namespace shimmer;

void
make_init_graph(infrastructure_graph& igraph)
{

    std::vector<vertex_descriptor> vds;

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( boost::add_vertex( { "station 0", 0, 0., 0.,  0., x}, igraph) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 0.,  0., x}, igraph) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 0.,  0., x}, igraph) );

    edge_properties ep0  = {edge_type::pipe, 0,    80000, 0.6, 0.000012};
    edge_properties ep1  = {edge_type::pipe, 1,    90000, 0.6, 0.000012};
    edge_properties ep2  = {edge_type::pipe, 2,   100000, 0.6, 0.000012};

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

    std::vector<double> ref_f = {9.037840034348319e-03,
                                9.167633023143456e-03,
                                1.113364421638781e-02};

    std::vector<double> ref_mu ={1.102400000000000e-05,
                                1.102400000000000e-05,
                                1.102400000000000e-05};

    bool pass= true;
    double Tm = 273.15 + 20.0; 

    size_t num_pipes = 3;
    size_t num_components = 17;

    vector_t flux (num_pipes);

    flux << 2.448496271836982e+01,
            2.171503728211028e+01,
            6.315037282326053e+00;

    infrastructure_graph graph;
    make_init_graph(graph);

    vector_t lambda(num_edges(graph));

    auto mu = viscosity<viscosity_type::Kukurugya>(Tm, graph); 

    size_t i = 0;
    auto edge_range = edges(graph);
    for(auto itor = edge_range.first; itor != edge_range.second; itor++,i++ )
    {
        auto pipe = graph[*itor];
        lambda(i) = friction_factor_average(pipe, Tm, flux(i), mu(i)); 
    }
    
    bool mu_pass = verify_test("Test viscosity", mu, ref_mu);
    bool f_pass = verify_test("Test friction factor", lambda, ref_f);
    
    return !(mu_pass && f_pass);
}