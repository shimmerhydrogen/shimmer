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

#include "MATLAB_GERG_functions.hpp"
#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/matlab_manip.h"
#include "../src/gas_law.h"
#include "verify_test.h"

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

    vds.push_back( boost::add_vertex({ "station 0", 0, 0., 0.,  0., x}, g));
    vds.push_back( boost::add_vertex({ "station 1", 1, 0., 0.,  0., x}, g));
    vds.push_back( boost::add_vertex({ "station 2", 2, 0., 0.,  0., x}, g));

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
    double Temp = 293.15;
    size_t num_pipes = 3;
    size_t num_nodes = 3;

    std::vector<double> ref_c2_nodes = {138267.2930151191,
                                        138578.7460692530,
                                        138546.3842273756};

    std::vector<double> ref_c2_pipes = {138422.1905008311,
                                        138406.1731482413,
                                        138562.5561693802};

    vector_t flux (num_pipes), flux_old(num_pipes);
    vector_t pressure_pipes (num_pipes), pressure_nodes(num_nodes);

    pressure_nodes << 5101325.0, 4977209.550248852, 4990077.876609823;     
    pressure_pipes << 5039522.018589198, 5045905.835430760, 4983646.482384357;

    vector_t RRp = make_RR(num_pipes); 
    vector_t RRn = make_RR(num_nodes); 

    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    gerg_params gerg_nodes = make_gerg(num_nodes); 
    gerg_params gerg_pipes = make_gerg(num_pipes); 

    auto  x_nodes = build_x_nodes(graph);
    Eigen::MatrixXd  x_pipes = inc.matrix_in().transpose() * x_nodes;

    auto eos_pipes = equation_of_state(Temp, pressure_pipes, x_pipes, gerg_pipes);
    auto eos_nodes = equation_of_state(Temp, pressure_nodes, x_nodes, gerg_nodes);

    vector_t c2_pipes = eos_pipes.Z.cwiseProduct(RRp) * Temp; 
    vector_t c2_nodes = eos_nodes.Z.cwiseProduct(RRn) * Temp; 

    std::cout << __FILE__ << std::endl; 
    bool c2n_pass = verify_test("Speed of sound at nodes", c2_nodes, ref_c2_nodes);
    bool c2p_pass = verify_test("Speed of sound in pipes", c2_pipes, ref_c2_pipes);
          
    return !(c2n_pass & c2p_pass);
}