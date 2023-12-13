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
#include "../src/conservation_matrices.h"
#include "../src/assemble.h"


using triple_t = std::array<double, 3>;

using namespace shimmer;


void
make_init_graph(infrastructure_graph& igraph)
{

    std::vector<vertex_descriptor> vds;

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
 
    boost::add_edge( vds[0], vds[1], ep0, igraph);
    boost::add_edge( vds[3], vds[1], ep1, igraph);
    boost::add_edge( vds[1], vds[2], ep2, igraph);
}


bool verify_test(const std::string & name, 
                 const sparse_matrix_t& mat,
                 const std::vector<triple_t>& ref )
{
    using itor_t = sparse_matrix_t::InnerIterator;
    bool pass = true;
    size_t count = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (itor_t it(mat,k); it; ++it, count++)
        { 
            std::cout << std::setprecision(16) << "(" << it.row() 
                      << " , " << it.col() << " , " << it.value() 
                      << " ) " << std::endl ;

            auto t = ref.at(count);
            auto e_val = std::abs((it.value() - t[2])/t[2]);
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
    std::vector<triple_t> ref_lhs ={{{0, 0, 9.621127501618740e-03},                                    
                                    {4, 0, 1}, 
                                    {1, 1, 1.350884841043611e-02},
                                    {4, 1,-0.503586391306371},
                                    {5, 1,-0.612626394184416},
                                    {6, 1, 1},
                                    {2 ,2, 2.474004214701962e-03},
                                    {6, 2,-1.341783903666971},
                                    {3, 3, 1.413716694115407e-03},
                                    {5, 3, 1},
                                    {0, 4, 1},
                                    {1, 4,-1},
                                    {4, 4, 6.763108109027953e+05}, 
                                    {1, 5,-1},
                                    {3, 5, 1},
                                    {5, 5, 4.558672924222292e+07},
                                    {1, 6, 1},
                                    {2, 6,-1}, 
                                    {6, 6, 1.173932795107726e+07}                                    
                                     }};


    double dt = 0.1;
    double c2 = 1000;

    infrastructure_graph graph;
    make_init_graph(graph);

    vector_t flux (num_edges(graph));
    vector_t pressure (num_vertices(graph)); 

    flux <<  -11, 13, -17; 
    pressure << 2000, 3000, 5000, 7000; 

    //incidence<double> inc(graph);
    assembler asmr(graph); 
    vector_t pm (num_edges(graph)); 
    average(pressure, asmr.incidence_in(), asmr.incidence_out(), pm);

    asmr.conservation_matrices(dt, c2, flux, pm, graph);
    asmr.assemble(graph);

    bool lhs_pass = verify_test("LHS matrix", asmr.LHS_matrix(), ref_lhs);

    return 0; 
}