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



template<typename GRAPH>
static void
make_init_graph(GRAPH& igraph)
{

 std::vector<typename GRAPH::vertex_descriptor> vds;

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



bool verify_test(const std::string & name, 
                 const vector_t<double>& vals,
                 const std::array<double, 3>& ref )
{
    using itor_t = Eigen::SparseMatrix<double>::InnerIterator;
    bool pass = true;
    size_t count = 0;
    for (int k = 0; k < vals.size(); ++k)
    {
        // std::cout << vals[k]  <<  std::endl ;

        auto e_val = (std::abs(vals[k] - ref.at(k))) /  ref.at(k);
        if(e_val > 1.e-12)
        {
            pass = false;
            break;  
        }
    }
    
    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << "  Test " << name << ".........." <<  passfail(pass) << std::endl;

    return pass;
}


int main()
{
    // Not realisic speed of sound. Intendeed only for test purposes.
    double c2 = 1; 
    double dt = 0.1; 

    std::array<double, 3> ref_pm = {2533.333333333333,5266.666666666666, 4083.333333333333};

    infrastructure_graph graph;
    make_init_graph(graph);

    Eigen::SparseMatrix<double> incidence_out = incidence_matrix_out<double>(graph);
    Eigen::SparseMatrix<double> incidence_in  = incidence_matrix_in<double>(graph);

    vector_t<double> pressure (num_vertices(graph)); 
    vector_t<double> pm (num_edges(graph)); 

    pressure << 2000, 3000, 5000, 7000; 

    average(pressure, incidence_in, incidence_out, pm);

    bool pass = verify_test(" mean pressure in pipes ", pm, ref_pm);

    return  !pass;

}