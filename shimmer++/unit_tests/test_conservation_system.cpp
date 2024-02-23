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
#include "verify_test.h"
#include "MATLAB_GERG_functions.hpp"
#include "../src/viscosity.h"


using triple_t = std::array<double, 3>;

using namespace shimmer;


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

    std::vector<triple_t> ref_lhs_mass =
                          {{{0 , 0 , 0.0009656491051934202}, 
                            {1 , 1 , 0.001020154052634392}, 
                            {2 , 2 , 0.001077080804942649}, 
                            {0 , 3 , 1}, 
                            {1 , 3 ,-1}, 
                            {0 , 4 , 1},
                            {2 , 4 ,-1},
                            {1 , 5 ,-1}, 
                            {2 , 5 , 1}}}; 
    std::vector<triple_t> ref_lhs_mom =
                          {{{3 , 0 , 1.0},
                            {4 , 0 , 1.0}, 
                            {3 , 1 ,-1.0}, 
                            {5 , 1 ,-1.0}, 
                            {4 , 2 ,-1.0},                             
                            {5 , 2 , 1.0}, 
                            {3 , 3 ,-11710.07601043195}, 
                            {4 , 4 ,-12014.55083082719}, 
                            {5 , 5 ,-6040.332315032241}}};

    std::vector<double> ref_rhs_mass = {4926.0899215508240,
                                        5077.5204934969630,
                                        5374.7170960654110};
    std::vector<double> ref_rhs_mom = {-1.626053247942385e+05,
                                       -1.496492958288440e+05,
                                       -2.527659740697931e+04};


    size_t num_pipes = 3;
    size_t num_nodes = 3;
    size_t num_edges_ext = num_pipes;
    size_t system_size = num_nodes + num_pipes;

    double dt = 180;
    double temperature = 293.15;

    vector_t flux (num_pipes), flux_old(num_pipes);
    vector_t pressure (num_nodes), pressure_old(num_nodes);
    vector_t c2_nodes (num_nodes), c2_pipes(num_pipes);  
 
    flux << 2.448496272217528e+01,2.171503727782473e+01, 6.315037277824726e+00; 
    c2_nodes << 138267.2930151191,138578.7460692530,138546.3842273756;
    c2_pipes << 138422.1905008311,138406.1731482413,138562.5561693802;
    pressure << 5101325.0, 4977209.550248852, 4990077.876609823;   
    pressure_old = pressure;
    flux_old = flux;
    //flux_ext = flux;

    bool unsteady = true;
    infrastructure_graph graph;
    make_init_graph(graph);

    auto mu = viscosity<viscosity_type::Kukurugya>(temperature, graph); 

    incidence inc(graph);
    linearized_fluid_solver lfs(unsteady, 0, dt, temperature,mu, inc, graph);

    vector_t pressure_pipes = average(pressure, inc);
    auto mass = lfs.continuity(pressure_old, c2_nodes);
    auto mom  = lfs.momentum(pressure, pressure_pipes, flux, flux_old, c2_pipes);

    sparse_matrix_t LHS_mass(num_nodes, system_size);
    sparse_matrix_t LHS_mom(system_size, system_size);

    LHS_mass.setFromTriplets(mass.first.begin(), mass.first.end());
    LHS_mom.setFromTriplets(mom.first.begin(), mom.first.end());

    bool lhs_mass_pass = verify_test("LHS continuity", LHS_mass , ref_lhs_mass);
    bool lhs_mom_pass  = verify_test("LHS momentum", LHS_mom, ref_lhs_mom);

    bool rhs_mass_pass = verify_test("rhs continuity", mass.second, ref_rhs_mass);
    bool rhs_mom_pass  = verify_test("rhs momentum", mom.second, ref_rhs_mom);

    return !(lhs_mass_pass && lhs_mom_pass && rhs_mass_pass && rhs_mom_pass); 
}