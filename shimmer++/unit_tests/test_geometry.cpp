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
#include "../src/geometry_properties.h"


static void
make_init_infrastructure(infrastructure_graph& igraph)
{

 std::vector<vertex_descriptor> vds;

    vds.push_back( boost::add_vertex( { "station 0", 0, 5000., -60, 0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 20 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 25 ,0. }, igraph) );

    edge_properties ep0  = {edge_type::pipe, 0,   5, 0.7, 0.012};
    edge_properties ep1  = {edge_type::pipe, 1,   9, 0.2, 0.012};
    edge_properties ep2  = {edge_type::pipe, 2,   7, 0.3, 0.012};
    edge_properties ep3  = {edge_type::pipe, 3,   2, 0.5, 0.012};

    /*           0                                *0  *1  *2  *3   
    //           |                              ----------------   
    //           |*0                           0|  1           
    //           1                             1| -1   1  -1   
    //         / |                             2|          1  -1
    //        /  |                             3|     -1       1   
    //     *2/   |*1                               
    //      /    |                             
    //   2 /_____3                              
             *3
    */

    boost::add_edge( 0, 1, ep0, igraph);
    boost::add_edge( 1, 3, ep1, igraph);
    boost::add_edge( 2, 1, ep2, igraph);
    boost::add_edge( 3, 2, ep3, igraph);
}



bool test(const std::array<double, 4>& ref)
{
    infrastructure_graph graph;
    make_init_infrastructure(graph);

    int i = 0;
    bool pass = true;
    auto v_range = vertices(graph);

    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
    {   
        //std::cout  << volume(*itor, igraph) << std::endl;

        if( std::abs(ref.at(i) - volume(*itor, graph)) > 1.e-12)
            pass = false;
    }

    return pass; 
}


int main(int argc, char **argv)
{

    std::array<double, 4> ref = {0.962112750161874, 1.350884841043611,
                                 0.443749962319558, 0.337721210260903} ; 

    bool pass = test(ref);

    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << __FILE__ << std::endl;
    std::cout << "  Test geometry ......" <<  passfail(pass) << std::endl;
    
    return !pass; 
}