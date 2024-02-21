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

matrix_t 
make_mass_fraction(size_t size)
{
    matrix_t mass_frac(size, 21);
    mass_frac.col(0).setConstant(1);

    return  mass_frac; 
}


gerg_params
make_gerg(size_t size)
{
    gerg_reducing_params_t reducing_parameters;
    reducing_parameters.Tr.resize(size,1);
    reducing_parameters.Dr.resize(size,1);
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

    size_t num_inlet = 1;    
    size_t num_pipes = 3;
    size_t num_nodes = 3;
    size_t num_bcnd = num_nodes;
    size_t system_size = num_nodes + num_pipes;

    double dt = 180;
    double temperature = 293.15;
    double tolerance = 1e-4;
    vector_t pressure(num_nodes), flux(num_pipes), L_rate(num_bcnd), flux_ext(num_pipes);
    
    flux_ext << 0.0, 30.80, 15.4;
    pressure << 5101325.0,
                4977209.550248852,
                4990077.876609823;
    flux <<     2.448496272217528e+01,
                2.171503727782473e+01,
                6.315037277824726e+00;
    L_rate <<   0.0, 0.0, 0.0;

    variable var(pressure, flux, L_rate);

    double pressure_in = 5101325.0;   
    vector_t inlet_nodes(num_inlet);
    inlet_nodes << 0; 


    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);
   
    matrix_t y_nodes = make_mass_fraction(num_nodes);
    matrix_t y_pipes = inc.matrix_in().transpose() * y_nodes;    

    vector_t area_pipes = area(graph);

    bool unsteady = true;

    gerg gerg_eos; 
    gerg_eos.compute_molar_mass(y_nodes, y_pipes);

    linearized_fluid_solver lfs(unsteady,tolerance, dt,temperature,inc, graph);
    lfs.run(area_pipes, inlet_nodes, pressure_in, flux_ext, var, var, &gerg_eos);

    vector_t sol(num_bcnd + num_pipes + num_nodes);
    sol.head(num_nodes) = var.pressure;
    sol.segment(num_nodes, num_pipes) = var.flux;
    sol.tail(num_bcnd) = var.L_rate;

    bool pass = verify_test("Test fluid-dynamic solver", sol, ref_sol); 

    return !pass;
}