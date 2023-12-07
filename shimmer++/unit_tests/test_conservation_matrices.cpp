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

 std::vector< typename GRAPH::vertex_descriptor> vds;

    vds.push_back( boost::add_vertex( { "station 0", 0, 0., 0., 100 }, igraph) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 0.,  30 }, igraph) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 0.,  60 }, igraph) );
    vds.push_back( boost::add_vertex( { "station 3", 3, 0., 0.,  80 }, igraph) );

    edge_properties ep0  = {edge_type::pipe, 0,   5, 0.7, 0.017};
    edge_properties ep1  = {edge_type::pipe, 1,   9, 0.2, 0.013};
    edge_properties ep2  = {edge_type::pipe, 2,   7, 0.3, 0.023};

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
                 const sparse_matrix_t<double>& mat,
                 const std::vector<triple_t>& ref )
{
    using itor_t = Eigen::SparseMatrix<double>::InnerIterator;
    bool pass = true;
    size_t count = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (itor_t it(mat,k); it; ++it, count++)
        {
            //std::cout << std::setprecision(16) << "(" << it.row() << " , " << it.col() << " , " <<  it.value()  << " )" <<  std::endl ;

            auto t = ref.at(count);
            auto e_val = std::abs(it.value() - t[2])/t[2];
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

    std::cout << "  Test " << name << " .........." <<  passfail(pass) << std::endl;

    return pass;
}

int main()
{
    std::vector<triple_t> ref_adp = {{{0,0,1}, {1,0, -0.503586391306371},
                                    {1, 1, -0.612626394184416},{3, 1, 1},
                                    {1, 2, 1}, {2, 2, -1.341783903666971}}};
    std::vector<triple_t> ref_resist = {{{0,0,6.763108109027953e+05}, 
                                    {1, 1, 4.558672924222292e+07}, 
                                    {2, 2, 1.173932795107726e+07}}};
    std::vector<triple_t> ref_phi = {{{0,0, 9.621127501618740e-03},
                                     {1, 1, 1.350884841043611e-02},
                                     {2 ,2, 2.474004214701962e-03},
                                     {3, 3, 1.413716694115407e-03}}};


    // Not realisic speed of sound. Intendeed only for test purposes.
    double c2 = 1000; 
    double dt = 0.1; 

    infrastructure_graph graph;
    make_init_graph(graph);

    vector_t<double> flux (num_edges(graph));
    vector_t<double> pressure (num_vertices(graph)); 

    flux <<  -11, 13, -17; 
    pressure << 2000, 3000, 5000, 7000; 

    sparse_matrix_t<double> incidence_out = incidence_matrix_out<double>(graph);
    sparse_matrix_t<double> incidence_in  = incidence_matrix_in<double>(graph);
    sparse_matrix_t<double> sADP(num_vertices(graph), num_edges(graph));
    sparse_matrix_t<double> sR(num_vertices(graph), num_edges(graph));
    sparse_matrix_t<double> sPHI(num_vertices(graph), num_vertices(graph));


    vector_t<double> pm (num_edges(graph)); 
    average(pressure, incidence_in, incidence_out, pm);

    adp_matrix(c2, graph, incidence_in, incidence_out, sADP);
    resistance_matrix(dt, c2, flux, pm, graph, sR);
    phi_matrix( dt,  c2, graph, sPHI);

    std::cout << __FILE__ << std::endl;

    bool adp_pass = verify_test("ADP matrix", sADP, ref_adp);
    bool res_pass = verify_test("R   matrix", sR  , ref_resist);
    bool phi_pass = verify_test("PHI matrix", sPHI, ref_phi);

    return !(adp_pass && res_pass && phi_pass) ;
}

