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
#include "verify_test.h"
#include "MATLAB_GERG_functions.hpp"

using triple_t = std::array<double, 3>;

using namespace shimmer;



vector_t
make_RR(size_t size)
{
    vector_t rrb(size);
    rrb.setConstant(518.2783563119373);
    return rrb; 
}


gerg_params
make_gerg(size_t size)
{
    gerg_reducing_params_t reducing_parameters;
    reducing_parameters.Tr.resize(size);
    reducing_parameters.Dr.resize(size);
    reducing_parameters.Tr.setConstant(1.905640000000000e+02);
    reducing_parameters.Dr.setConstant(1.013934271900000e+01);

    gerg_pseudo_critical_pt_t psc_point;
    psc_point.Tcx.resize(size, 1);
    psc_point.Dcx.resize(size, 1);
    psc_point.Vcx.resize(size, 1);
    psc_point.Tcx.setConstant(1.905640000000000e+02);
    psc_point.Dcx.setConstant(1.013934271900000e+01);
    psc_point.Vcx.setConstant(9.862572236818776e-02);
    
    gerg_thermo_params_t parameters;
    parameters.Type = gerg_thermo_params_t::Types::Gas_phase;
    
    gerg_params gerg(reducing_parameters, psc_point, parameters); 

    return gerg;
}



void
make_init_graph(infrastructure_graph& g)
{

    std::vector<vertex_descriptor> vds;

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( boost::add_vertex( { "station 0", 0, 0., 0.,  0., x}, g) );
    vds.push_back( boost::add_vertex( { "station 1", 1, 0., 0.,  0., x}, g) );
    vds.push_back( boost::add_vertex( { "station 2", 2, 0., 0.,  0., x}, g) );

    edge_properties ep0  = {edge_type::pipe, 0,    80000, 0.6, 0.000012};
    edge_properties ep1  = {edge_type::pipe, 1,    90000, 0.6, 0.000012};
    edge_properties ep2  = {edge_type::pipe, 2,   100000, 0.6, 0.000012};

    /*                                            
    //           0                        *0  *1  *2              
    //         / |                     --------------         
    //        /  |                     0|  1   1                        
    //     *1/   |*0                   1| -1      -1          
    //      /____|                     2|     -1   1        
    //    1   *2   2                                         
    */

    boost::add_edge( vds[0], vds[1], ep0, g);
    boost::add_edge( vds[0], vds[2], ep1, g);
    boost::add_edge( vds[2], vds[1], ep2, g);
}



int main()
{

    std::vector<triple_t> ref_lhs =
                          {{{0 , 0 , 0.0009656491051934202}, 
                            {3 , 0 , 10078534.55024885},
                            {4 , 0 , 10091402.87660982}, 
                            {1 , 1 , 0.001020154052634392}, 
                            {3 , 1 ,-10078534.55024885}, 
                            {5 , 1 ,-9967287.426858675}, 
                            {2 , 2 , 0.001077080804942649}, 
                            {4 , 2 ,-10091402.87660982},                             
                            {5 , 2 , 9967287.426858675}, 
                            {0 , 3 , 1}, 
                            {1 , 3 ,-1}, 
                            {3 , 3 ,-118020405657.1786}, 
                            {0 , 4 , 1},
                            {2 , 4 ,-1},
                            {4 , 4 ,-121243672815.3844}, 
                            {1 , 5 ,-1}, 
                            {2 , 5 , 1}, 
                            {5 , 5 ,-60205728337.66899}}};

    std::vector<double> ref_rhs = { 4926.0899215508240,
                                    5077.5204934969630,
                                    5374.7170960654110,
                                    -1638823383993.169,
                                    -1510171334409.830,
                                    -251939111528.3534};


    size_t num_pipes = 3;
    size_t num_vertices = 3;
    size_t num_edges_ext = num_pipes;

    double dt = 180;
    double temperature = 293.15;

    vector_t flux (num_pipes), flux_old(num_pipes);
    vector_t pressure (num_vertices), pressure_old(num_vertices);
    vector_t c2_vertices (num_vertices), c2_edges(num_pipes); 

    flux << 2.448496272217528e+01,
            2.171503727782473e+01,
            6.315037277824726e+00; 

    flux_old = flux;
    //flux_ext = flux;

    pressure << 5101325.0, 4977209.550248852, 4990077.876609823;
    
    pressure_old = pressure;

    c2_vertices << 138267.2930151191,
                138578.7460692530,
                138546.3842273756;

    c2_edges << 138422.1905008311,
                138406.1731482413,
                138562.5561693802;
    
    gerg_params gerg_nodes = make_gerg(num_vertices); 
    gerg_params gerg_pipes = make_gerg(num_pipes); 

    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    auto x_nodes = build_x_nodes(graph);
    auto x_pipes = inc.matrix_in().transpose() * x_nodes;
    auto RR_nodes = make_RR(num_vertices);
    auto RR_pipes = make_RR(num_pipes);

    auto system_mass = continuity(dt, temperature, pressure, pressure_old,
                                  inc, graph, x_nodes, RR_nodes, gerg_nodes);
    auto system_mom  = momentum(dt, temperature, flux, flux_old, pressure, 
                                  inc, graph, x_pipes, RR_pipes, gerg_pipes);

    auto [LHS, rhs] = assemble(system_mass, system_mom, graph);

    /*
    std::cout << "LHS: " <<  std::endl ;
    size_t count = 0;
    for (int k = 0; k < LHS.outerSize(); ++k)
    {
        for (itor_t it(LHS,k); it; ++it, count++)
        { 
            std::cout << std::setprecision(16) << "{" << it.row() 
                      << " , " << it.col() << " , " << it.value() 
                      << "} " << std::endl ;
        }
    }

    std::cout << "rhs: " <<  std::endl ;
    for (int k = 0; k < rhs.size(); ++k)
        std::cout <<  rhs(k)<< std::endl ;
    std::cout <<  std::endl ;
    */  

    bool lhs_pass = verify_test("LHS matrix", LHS, ref_lhs);
    bool rhs_pass = verify_test("rhs vector", rhs, ref_rhs);

    return !(lhs_pass && rhs_pass); 
}