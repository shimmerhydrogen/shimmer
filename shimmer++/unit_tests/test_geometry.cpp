/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "solver/geometry_properties.h"
#include "verify_test.h"


using namespace shimmer;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


static void
make_init_infrastructure(infrastructure_graph& igraph)
{

    std::vector<vertex_descriptor> vds;


    auto add_vertex = [&](vertex_properties&& vp) 
    {
        auto v = boost::add_vertex(igraph);
        igraph[v] = std::move(vp);
        return v;
    };

    vds.push_back( add_vertex( vertex_properties( "station 0", 0, 5000., -60, 0. ) ) );
    vds.push_back( add_vertex( vertex_properties(  "station 1", 1, 0., 20 ,0. )) );
    vds.push_back( add_vertex( vertex_properties(  "station 2", 2, 0., 25 ,0. )) );
    vds.push_back( add_vertex( vertex_properties(  "station 3", 2, 0., 25 ,0. )) );

    edge_properties ep0  = {pipe_type::PIPE, 0,   5, 0.7, 0.012};
    edge_properties ep1  = {pipe_type::PIPE, 1,   9, 0.2, 0.012};
    edge_properties ep2  = {pipe_type::PIPE, 2,   7, 0.3, 0.012};
    edge_properties ep3  = {pipe_type::PIPE, 3,   2, 0.5, 0.012};

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

    boost::add_edge( vds[0], vds[1], ep0, igraph);
    boost::add_edge( vds[1], vds[3], ep1, igraph);
    boost::add_edge( vds[2], vds[1], ep2, igraph);
    boost::add_edge( vds[3], vds[2], ep3, igraph);
}


int main(int argc, char **argv)
{

    std::vector<double> ref = {0.962112750161874, 1.350884841043611,
                                 0.443749962319558, 0.337721210260903} ; 

    infrastructure_graph graph;
    make_init_infrastructure(graph);

    vector_t vols (num_vertices(graph));

    size_t i = 0;
    auto v_range = vertices(graph);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        vols(i) = volume(*itor, graph);

    std::cout << __FILE__ << std::endl;
    bool pass = verify_test("geometry", vols, ref);    
    return !pass; 
}