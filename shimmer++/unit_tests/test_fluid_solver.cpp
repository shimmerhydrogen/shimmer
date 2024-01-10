/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
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
#include "../src/fluid_solver.h"

using triple_t = std::array<double, 3>;

using namespace shimmer;


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


std::pair<vector_t, vector_t>
make_rr_mm(size_t size)
{
    vector_t rrb(size);
    rrb.setConstant(518.2783563119373);

    vector_t molar_mass(size);
    molar_mass.setConstant(16.04246);

    return std::make_pair(rrb, molar_mass); 
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


int main()
{
     std::vector<double> ref_sol = {5101325,
                                    4977209.550249534,
                                    4990077.876609622,
                                    24.48496271836996,
                                    21.71503728211033,
                                    6.315037282326063,
                                    -46.20000000048003,
                                    30.8,
                                    15.4};

    size_t num_pipes = 3;
    size_t num_nodes = 3;
    size_t num_bcnd = num_nodes;
    size_t system_size = num_nodes + num_pipes;

    double dt = 180;
    double temperature = 293.15;
    double tolerance = 1e-4;
    vector_t sol (num_pipes + num_nodes + num_bcnd), flux_ext(num_pipes);
    
    flux_ext << 0.0, 30.80, 15.4;

    sol <<  5101325.0,
            4977209.550248852,
            4990077.876609823,
            2.448496272217528e+01,
            2.171503727782473e+01,
            6.315037277824726e+00,
            0.0, 0.0, 0.0;

    double pressure_in = 5101325.0; 

    gerg_params gerg_nodes = make_gerg(num_nodes); 
    gerg_params gerg_pipes = make_gerg(num_pipes); 

    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    auto x_nodes = build_x_nodes(graph);
    auto x_pipes = inc.matrix_in().transpose() * x_nodes;
    auto [rr_nodes, mm_nodes] = make_rr_mm(num_nodes);
    auto [rr_pipes, mm_pipes] = make_rr_mm(num_pipes);


    linearized_fluid_solver( tolerance, dt,temperature,
                        inc, graph, pressure_in,
                        flux_ext, rr_nodes,
                        rr_pipes, mm_pipes, gerg_nodes, gerg_pipes,
                        sol);

    bool pass = verify_test("Test fluid-dynamic solver", sol, ref_sol); 

    return !pass;
}