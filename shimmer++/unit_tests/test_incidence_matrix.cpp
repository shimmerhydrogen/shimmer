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

#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "verify_test.h"


using namespace shimmer;


static void
make_init_graph(infrastructure_graph& igraph)
{


    /*          0                                  0  1  2  3  4  5  6  7  8  9  10  11   
    //           |                              ----------------------------------------   
    //           |*0                           0|  1           
    //           1                             1| -1  1     1   
    //         / |                             2|        1 -1     1
    //        /  |                             3|                -1  1        1      -1  
    //     *3/   |*1                           4|    -1 -1  1              1          1    
    //      /    |                             5|             -1        1    -1   1
    //   2 /_____4_____7                       6|                   -1 -1 
    //    |  *2 /| *8  /                       7|                         -1     -1
    //    |    / |    /
    //    |  /   |*4 / 
    //  *5| /*11 |  /*10 
    //    |/_____| /
    //   3|  *9 / 5
    //    |    /
    //  *6|   /*7 
    //    |  /  
    //    | /      
          6 
    */

    std::vector<vertex_descriptor> vds;

    auto add_vertex = [&](vertex_properties&& vp) 
    {
        auto v = boost::add_vertex(igraph);
        igraph[v] = std::move(vp);
        return v;
    };


    vds.push_back(add_vertex( vertex_properties("station 0", 0, 5000., -60, 0.)) );
    vds.push_back(add_vertex( vertex_properties( "station 1", 1, 0., 20 ,0. )) );
    vds.push_back(add_vertex( vertex_properties( "station 2", 2, 0., 25 ,0. )) );
    vds.push_back(add_vertex( vertex_properties( "station 3", 3, 0., 30 ,0. )) );
    vds.push_back(add_vertex( vertex_properties( "station 4", 4, 0., 35 ,0. )) );
    vds.push_back(add_vertex( vertex_properties( "station 5", 5, 0., 40 ,0. )) );
    vds.push_back(add_vertex( vertex_properties( "station 6", 6, 0., 45 ,0. )) );
    vds.push_back(add_vertex( vertex_properties( "station 7", 7, 0., 50 ,0. )) );

    using eprop_t = edge_properties;

    edge_properties ep0  = {pipe_type::PIPE, 0,   80, 0.6, 0.012};
    edge_properties ep1  = {pipe_type::PIPE, 1,   90, 0.6, 0.012};
    edge_properties ep2  = {pipe_type::PIPE, 2,  100, 0.6, 0.012};
    edge_properties ep3  = {pipe_type::PIPE, 3,  110, 0.6, 0.012};
    edge_properties ep4  = {pipe_type::PIPE, 4,   80, 0.6, 0.012};
    edge_properties ep5  = {pipe_type::PIPE, 5,   80, 0.6, 0.012};
    edge_properties ep6  = {pipe_type::PIPE, 6,   80, 0.6, 0.012};
    edge_properties ep7  = {pipe_type::PIPE, 7,   80, 0.6, 0.012};
    edge_properties ep8  = {pipe_type::PIPE, 8,   80, 0.6, 0.012};
    edge_properties ep9  = {pipe_type::PIPE, 9,   80, 0.6, 0.012};
    edge_properties ep10 = {pipe_type::PIPE,10,   80, 0.6, 0.012};
    edge_properties ep11 = {pipe_type::PIPE,11,   80, 0.6, 0.012};


    boost::add_edge( vds[0], vds[ 1], ep0, igraph);
    boost::add_edge( vds[1], vds[ 4], ep1, igraph);
    boost::add_edge( vds[2], vds[ 4], ep2, igraph);
    boost::add_edge( vds[1], vds[ 2], ep3, igraph);
    boost::add_edge( vds[4], vds[ 5], ep4, igraph);
    boost::add_edge( vds[2], vds[ 3], ep5, igraph);
    boost::add_edge( vds[3], vds[ 6], ep6, igraph);
    boost::add_edge( vds[5], vds[ 6], ep7, igraph);
    boost::add_edge( vds[4], vds[ 7], ep8, igraph);
    boost::add_edge( vds[3], vds[ 5], ep9, igraph);
    boost::add_edge( vds[5], vds[ 7], ep10, igraph);
    boost::add_edge( vds[4], vds[ 3], ep11, igraph);
}


int main(int argc, char **argv)
{
    using triple_t = std::array<double, 3>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>; 
    
    std::vector<triple_t> ref = {{{0,0,1},{1,0,-1},{1,1,1},{4,1,-1},{2,2,1},{4,2,-1},
                                 {1,3,1},{2,3,-1},{4,4,1},{5,4,-1},{2,5,1},{3,5,-1},
                                 {3,6,1},{6,6,-1},{5,7,1},{6,7,-1},{4,8,1},{7,8,-1},
                                 {3,9,1},{5,9,-1},{5,10,1},{7,10,-1},{3,11,-1},{4,11,1}}};

    infrastructure_graph igraph;
    make_init_graph(igraph);

    incidence inc(igraph);
    const sparse_matrix_t& mat  = inc.matrix();

    std::cout << __FILE__ << std::endl;
    bool pass = verify_test("incidence matrix", mat, ref);
    
    return !pass; 
}

