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


#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"

static void
make_init_infrastructure(infrastructure_graph& igraph)
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

    vds.push_back( boost::add_vertex( { "station 0", 0, 5000., -60, 0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 20 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 25 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 3", 3, 0., 30 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 4", 4, 0., 35 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 5", 5, 0., 40 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 6", 6, 0., 45 ,0. }, igraph) );
    vds.push_back( boost::add_vertex( { "station 7", 7, 0., 50 ,0. }, igraph) );


    using eprop_t = edge_properties;

    edge_properties ep0  = {edge_type::pipe, 0,   80, 0.6, 0.012};
    edge_properties ep1  = {edge_type::pipe, 1,   90, 0.6, 0.012};
    edge_properties ep2  = {edge_type::pipe, 2,  100, 0.6, 0.012};
    edge_properties ep3  = {edge_type::pipe, 3,  110, 0.6, 0.012};
    edge_properties ep4  = {edge_type::pipe, 4,   80, 0.6, 0.012};
    edge_properties ep5  = {edge_type::pipe, 5,   80, 0.6, 0.012};
    edge_properties ep6  = {edge_type::pipe, 6,   80, 0.6, 0.012};
    edge_properties ep7  = {edge_type::pipe, 7,   80, 0.6, 0.012};
    edge_properties ep8  = {edge_type::pipe, 8,   80, 0.6, 0.012};
    edge_properties ep9  = {edge_type::pipe, 9,   80, 0.6, 0.012};
    edge_properties ep10 = {edge_type::pipe,10,   80, 0.6, 0.012};
    edge_properties ep11 = {edge_type::pipe,11,   80, 0.6, 0.012};


    boost::add_edge( 0, 1, ep0, igraph);
    boost::add_edge( 1, 4, ep1, igraph);
    boost::add_edge( 2, 4, ep2, igraph);
    boost::add_edge( 1, 2, ep3, igraph);
    boost::add_edge( 4, 5, ep4, igraph);
    boost::add_edge( 2, 3, ep5, igraph);
    boost::add_edge( 3, 6, ep6, igraph);
    boost::add_edge( 5, 6, ep7, igraph);
    boost::add_edge( 4, 7, ep8, igraph);
    boost::add_edge( 3, 5, ep9, igraph);
    boost::add_edge( 5, 7, ep10, igraph);
    boost::add_edge( 4, 3, ep11, igraph);


}

int main(int argc, char **argv)
{
    using triple_t = std::array<int, 3>;
    std::array<triple_t, 24> ref = {{{0,0,1},{1,0,-1},{1,1,1},{4,1,-1},{2,2,1},{4,2,-1},
                                 {1,3,1},{2,3,-1},{4,4,1},{5,4,-1},{2,5,1},{3,5,-1},
                                 {3,6,1},{6,6,-1},{5,7,1},{6,7,-1},{4,8,1},{7,8,-1},
                                 {3,9,1},{5,9,-1},{5,10,1},{7,10,-1},{3,11,-1},{4,11,1}}};


    infrastructure_graph igraph;
    make_init_infrastructure(igraph);
    Eigen::SparseMatrix<int> mat = incidence_matrix<int>(igraph);

    /*
    std::cout << "Incidence matrix: \n"; 
    for (int k=0; k<mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat,k); it; ++it)
            std::cout << "{"<< it.row() << ","<< it.col() <<"," << it.value()<< "},";
    std::cout << std::endl; 
    */

    bool pass = true;
    size_t count = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat,k); it; ++it, count++)
        {
            auto t = ref.at(count);
            if((it.row() != t[0])  || (it.col() != t[1]) || (it.value() != t[2]))
            {
                pass = false;
                break;  
            }
        }
    }

    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << __FILE__ << std::endl;
    std::cout << "  Test incidence matrix ......" <<  passfail(pass) << std::endl;
    
    return pass; 
}

