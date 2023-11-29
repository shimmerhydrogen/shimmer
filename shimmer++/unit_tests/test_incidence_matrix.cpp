/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
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

    edge_properties ep0  = {edge_type::pipe, 0, 0, 1,  80, 0.6, 0.012};
    edge_properties ep1  = {edge_type::pipe, 1, 1, 4,  90, 0.6, 0.012};
    edge_properties ep2  = {edge_type::pipe, 2, 2, 4, 100, 0.6, 0.012};
    edge_properties ep3  = {edge_type::pipe, 3, 1, 2, 110, 0.6, 0.012};
    edge_properties ep4  = {edge_type::pipe, 4, 4, 5,  80, 0.6, 0.012};
    edge_properties ep5  = {edge_type::pipe, 5, 2, 3,  80, 0.6, 0.012};
    edge_properties ep6  = {edge_type::pipe, 6, 3, 6,  80, 0.6, 0.012};
    edge_properties ep7  = {edge_type::pipe, 7, 5, 6,  80, 0.6, 0.012};
    edge_properties ep8  = {edge_type::pipe, 8, 4, 7,  80, 0.6, 0.012};
    edge_properties ep9  = {edge_type::pipe, 9, 3, 5,  80, 0.6, 0.012};
    edge_properties ep10 = {edge_type::pipe,10, 5, 7,  80, 0.6, 0.012};
    edge_properties ep11 = {edge_type::pipe,11, 4, 3,  80, 0.6, 0.012};


    boost::add_edge( vds[ep0.from], vds[ep0.to], ep0, igraph);
    boost::add_edge( vds[ep1.from], vds[ep1.to], ep1, igraph);
    boost::add_edge( vds[ep2.from], vds[ep2.to], ep2, igraph);
    boost::add_edge( vds[ep3.from], vds[ep3.to], ep3, igraph);
    boost::add_edge( vds[ep4.from], vds[ep4.to], ep4, igraph);
    boost::add_edge( vds[ep5.from], vds[ep5.to], ep5, igraph);
    boost::add_edge( vds[ep6.from], vds[ep6.to], ep6, igraph);
    boost::add_edge( vds[ep7.from], vds[ep7.to], ep7, igraph);
    boost::add_edge( vds[ep8.from], vds[ep8.to], ep8, igraph);
    boost::add_edge( vds[ep9.from], vds[ep9.to], ep9, igraph);
    boost::add_edge( vds[ep10.from], vds[ep10.to], ep10, igraph);
    boost::add_edge( vds[ep11.from], vds[ep11.to], ep11, igraph);


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

