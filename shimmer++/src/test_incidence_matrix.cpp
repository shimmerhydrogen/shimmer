/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#include <iostream>
#include <fstream>
#include <string>


#include "infrastructure_graph.h"
#include "incidence_matrix.h"

static void
make_init_infrastructure(infrastructure_graph& igraph)
{
    std::vector<vertex_descriptor> vds;

    vds.push_back( boost::add_vertex( { "station0", 0, 5000., -60, 0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station1", 1, 0., 20 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station2", 2, 0., 25 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station3", 3, 0., 30 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station4", 4, 0., 35 ,0. }, igraph) );

    edge_properties ep0 = {edge_type::pipe, 0, 0, 1,  80, 0.6, 0.012};
    edge_properties ep1 = {edge_type::pipe, 1, 1, 2,  90, 0.6, 0.012};
    edge_properties ep2 = {edge_type::pipe, 2, 2, 3, 100, 0.6, 0.012};
    edge_properties ep3 = {edge_type::pipe, 3, 1, 4, 110, 0.6, 0.012};

    boost::add_edge( vds[0], vds[1], ep0, igraph);
    boost::add_edge( vds[1], vds[2], ep1, igraph);
    boost::add_edge( vds[2], vds[3], ep2, igraph);
    boost::add_edge( vds[1], vds[4], ep3, igraph);
}

int main(int argc, char **argv)
{
    infrastructure_graph igraph;

    make_init_infrastructure(igraph);

    Eigen::SparseMatrix<double> mat = incidence_matrix(igraph);

    std::cout << "Incidence matrix: \n"; 
    for (int k=0; k<mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
            std::cout << "("<< it.row() << ";"<< it.col() <<") = " << it.value() << std::endl;

    return 0;
}

