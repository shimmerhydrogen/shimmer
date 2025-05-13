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
#include "solver/pipe_calculator.h"
#include "solver/viscosity.h"
#include "verify_test.h"

using namespace shimmer;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


static void
make_init_graph(infrastructure_graph& igraph)
{

    std::vector<vertex_descriptor> vds;

    auto add_vertex = [&](vertex_properties&& vp, const vector_t& x_in) 
    {
        vp.gas_mixture = x_in;
        auto v = boost::add_vertex(igraph);
        igraph[v] = std::move(vp);
        return v;
    };

    vector_t x(21);
    x.setConstant(1.);


    vds.push_back( add_vertex(vertex_properties( "station 0", 0, 5000.,-60,0.), x));
    vds.push_back( add_vertex(vertex_properties( "station 1", 1, 0., 20 ,0.), x));
    vds.push_back( add_vertex(vertex_properties( "station 2", 2, 0., 25 ,0.), x));
    vds.push_back( add_vertex(vertex_properties( "station 3", 3, 0., 25 ,0.), x));

    edge_properties ep0  = {pipe_type::PIPE, 0,   5, 0.7, 0.00001};
    edge_properties ep1  = {pipe_type::PIPE, 1,   9, 0.2, 0.00003};
    edge_properties ep2  = {pipe_type::PIPE, 2,   7, 0.3, 0.00005};
    edge_properties ep3  = {pipe_type::PIPE, 3,   2, 0.5, 0.00007};
    edge_properties ep4  = {pipe_type::PIPE, 4,   1, 0.1, 0.00013};


    /*          _ 0                                *0  *1  *2  *3  *4  
    //         /  |                              ---------------------   
    //     *4 /   |*0                           0|  1              -1
    //       /    1                             1| -1   1  -1   
    //      /   / |                             2|          1  -1   1
    //     /   /  |                             3|     -1       1   
    //    |   /   |*1                               
    //    |  /*2  |                             
    //    2 /_____3                              
    //       *3
    */

    boost::add_edge( vds[0], vds[1], ep0, igraph);
    boost::add_edge( vds[1], vds[3], ep1, igraph);
    boost::add_edge( vds[2], vds[1], ep2, igraph);
    boost::add_edge( vds[3], vds[2], ep3, igraph);
    boost::add_edge( vds[2], vds[0], ep4, igraph);
}


int main(int argc, char **argv)
{
    double dt = 0.1;
    double c2 = 1.0;
    double p = 1.0;
    double temperature = 293.15;


    std::vector<double> ref_inertia  = {2.598448050479924e+02,
                                        5.729577951308232e+03,
                                        1.980594847365808e+03,
                                        2.037183271576260e+02,
                                        2.546479089470325e+03}; 
    std::vector<double> ref_friction = {7.586693363446199e-01,
                                        5.640811629945327e+02,
                                        5.823770872835033e+01,
                                        1.291280275074440e+00,
                                        3.372739194658072e+03}; 

    infrastructure_graph graph;
    make_init_graph(graph);  

    vector_t ri (num_edges(graph));
    vector_t rf (num_edges(graph));

    vector_t flux (num_edges(graph));
    flux << 19.0, 23.0, 37.0, 51.0, 17.0;

    auto mu = viscosity<viscosity_type::Kukurugya>(temperature, graph); 

    int i = 0;
    auto e_range = edges(graph);
    for(auto itor = e_range.first; itor != e_range.second; itor++, i++)
    {   
        auto pipe = graph[*itor];
        ri(i) = inertia_resistance(pipe, dt, p);
        rf(i) = friction_resistance(pipe, temperature, mu(i), c2, flux(i));
    }

    std::cout << __FILE__ << std::endl; 
    bool ipass = verify_test("resistance inertia", ri, ref_inertia);
    bool fpass = verify_test("resistance friction", rf, ref_friction);

    return !(ipass && fpass); 
}