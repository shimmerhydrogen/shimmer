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
#include "solver/conservation_matrices.h"
#include "verify_test.h"


using namespace shimmer;


static void
make_init_graph(infrastructure_graph& igraph)
{

    std::vector<vertex_descriptor> vds;

    auto add_vertex = [&](vertex_properties&& vp) 
    {
        auto v = boost::add_vertex(igraph);
        igraph[v] = std::move(vp);
        return v;
    };

    vds.push_back( add_vertex( vertex_properties( "station 0", 0, 0., 0., 100 )) );
    vds.push_back( add_vertex( vertex_properties( "station 1", 1, 0., 0.,  30 )) );
    vds.push_back( add_vertex( vertex_properties( "station 2", 2, 0., 0.,  60 )) );
    vds.push_back( add_vertex( vertex_properties( "station 3", 3, 0., 0.,  80 )) );

    edge_properties ep0  = {pipe_type::PIPE, 0,   5, 0.7, 0.012};
    edge_properties ep1  = {pipe_type::PIPE, 1,   9, 0.2, 0.012};
    edge_properties ep2  = {pipe_type::PIPE, 2,   7, 0.3, 0.012};

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