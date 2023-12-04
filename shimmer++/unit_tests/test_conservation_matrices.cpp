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
#include "../src/conservation_matrices.cpp"


using triple_t = std::array<double, 3>;


template<typename GRAPH>
static void
make_init_graph(GRAPH& igraph)
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

    boost::add_edge( 0, 1, ep0, igraph);
    boost::add_edge( 3, 1, ep1, igraph);
    boost::add_edge( 1, 2, ep2, igraph);
}


void verify_test(const std::string & name, 
                 const Eigen::SparseMatrix<double>& mat,
                 const std::array<triple_t, 6>& ref )
{
    using itor_t = Eigen::SparseMatrix<double>::InnerIterator;
    bool pass = true;
    size_t count = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (itor_t it(mat,k); it; ++it, count++)
        {
            std::cout << "(" << it.row() << " , " << it.col() << " , " <<  it.value()  << " )" <<  std::endl ;

            auto t = ref.at(count);
            auto e_val = std::abs(it.value() - t[2]);
            if((it.row() != t[0])  || (it.col() != t[1]) || (e_val > 1.e-12))
            {
                pass = false;
                break;  
            }
        }
    }
    
    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << "  Test " << name << "  matrix .........." <<  passfail(pass) << std::endl;

    return;
}

int main()
{
    // Not realisic speed of sound. Intendeed only for test purposes.
    double c2 = 1000; 

    std::array<triple_t, 6> ref_adp = {{{0,0,1}, {1,0, -0.503586391306371},
                                    {1, 1, -0.612626394184416},{3, 1, 1},
                                    {1, 2, 1}, {2, 2, -1.341783903666971}}};
    //std::array<triple_t, 6> ref_resist= {{{0,0,10}, {1,0,10},
    //                                {1, 1, 10}, {3, 1, 10},
    //                                {1, 2, 10}, {2, 2, 10}}};

    undirected_graph graph;
    make_init_graph(graph);

    Eigen::SparseMatrix<double> incidence_out = incidence_matrix_out<double>(graph);
    Eigen::SparseMatrix<double> incidence_in  = incidence_matrix_in<double>(graph);
    Eigen::SparseMatrix<double> sADP(num_vertices(graph), num_edges(graph));
    Eigen::SparseMatrix<double> sR(num_vertices(graph), num_edges(graph));

    adp_matrix(c2, graph, incidence_in, incidence_out, sADP);
    //resistance_matrix(c2, graph, incidence_in, incidence_out, sR);

    std::cout << __FILE__ << std::endl;

    verify_test("ADP", sADP, ref_adp);
    //verify_test("R", sR, ref_resist);


    return 0;
}

