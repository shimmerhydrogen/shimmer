/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "MATLAB_GERG_functions.hpp"
#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "solver/matlab_manip.h"
#include "solver/gas_law.h"
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
    reducing_parameters.Tr.resize(size, 1);
    reducing_parameters.Dr.resize(size, 1);
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

    auto add_vertex = [&](vertex_properties&& vp, const vector_t& x_in)
    {
        vp.gas_mixture = x_in;
        auto v = boost::add_vertex(g);
        g[v] = std::move(vp);
        return v;
    };

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( add_vertex( vertex_properties("station 0", 0, 0., 0.,  0.), x));
    vds.push_back( add_vertex( vertex_properties("station 1", 1, 0., 0.,  0.), x));
    vds.push_back( add_vertex( vertex_properties("station 2", 2, 0., 0.,  0.), x));

    edge_properties ep0  = {pipe_type::PIPE, 0,    80000, 0.6, 0.000012};
    edge_properties ep1  = {pipe_type::PIPE, 1,    90000, 0.6, 0.000012};
    edge_properties ep2  = {pipe_type::PIPE, 2,   100000, 0.6, 0.000012};

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


bool
gerg_matlab()
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
    vector_t temperature_nodes(num_nodes), temperature_pipes(num_pipes);

    pressure_nodes << 5101325.0, 4977209.550248852, 4990077.876609823;     
    pressure_pipes << 5039522.018589198, 5045905.835430760, 4983646.482384357;
    
    temperature_nodes.setConstant(Temp);
    temperature_pipes.setConstant(Temp);

    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    gerg_params gerg_nodes = make_gerg(num_nodes); 
    gerg_params gerg_pipes = make_gerg(num_pipes); 

    Eigen::MatrixXd  x_nodes = build_x_nodes(graph);
    Eigen::MatrixXd  x_pipes = inc.matrix_in().transpose() * x_nodes;

    gerg gerg_eos; 

    gerg_eos.compute_molar_mass(x_nodes, x_pipes);

    vector_t RRp = gerg_eos.Rgas_pipes(); 
    vector_t RRn = gerg_eos.Rgas_nodes(); 
    
    auto Z_nodes = gerg_eos.compute(Temp, pressure_nodes, x_nodes, gerg_nodes);
    auto Z_pipes = gerg_eos.compute(Temp, pressure_pipes, x_pipes, gerg_pipes);

    vector_t c2_nodes = Z_nodes.array() * RRn.array() * Temp; 
    vector_t c2_pipes = Z_pipes.array() * RRp.array() * Temp; 

    std::cout << __FILE__ << std::endl; 
    bool c2n_pass = verify_test("Speed of sound at nodes", c2_nodes, ref_c2_nodes);
    bool c2p_pass = verify_test("Speed of sound in pipes", c2_pipes, ref_c2_pipes);
          
    return !(c2n_pass & c2p_pass);
}

bool
RR_matlab()
{
    size_t num_pipes = 3;
    size_t num_nodes = 3;

    std::vector<double> ref_RRp = {518.2783563119373,
                                    518.2783563119373,
                                    518.2783563119373};

    std::vector<double> ref_RRn = {518.2783563119373,
                                    518.2783563119373,
                                    518.2783563119373};

    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    Eigen::MatrixXd  x_nodes = build_x_nodes(graph);
    Eigen::MatrixXd  x_pipes = inc.matrix_in().transpose() * x_nodes;

    gerg gerg_eos; 

    gerg_eos.compute_molar_mass(x_nodes, x_pipes);

    vector_t RRp = gerg_eos.Rgas_pipes();
    vector_t RRn = gerg_eos.Rgas_nodes(); 

    std::cout << __FILE__ << std::endl; 
    bool RRp_pass = verify_test("Rgas in pipes", RRp, ref_RRp);
    bool RRn_pass = verify_test("Rgas at nodes", RRn, ref_RRn);

    return !(RRp_pass & RRn_pass);
}


bool
gerg_aga8code()
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
    vector_t temperature_nodes(num_nodes), temperature_pipes(num_pipes);

    pressure_nodes << 5101325.0, 4977209.550248852, 4990077.876609823;     
    pressure_pipes << 5039522.018589198, 5045905.835430760, 4983646.482384357;

    temperature_pipes.setConstant(Temp);
    temperature_nodes.setConstant(Temp);

    
    infrastructure_graph graph;
    make_init_graph(graph);

    incidence inc(graph);

    matrix_t  x_nodes = build_x_nodes(graph);
    matrix_t  x_pipes = inc.matrix_in().transpose() * x_nodes;

    gerg_aga gerg_eos; 

    gerg_eos.compute_molar_mass(x_nodes, x_pipes);

    vector_t RRp = gerg_eos.Rgas_pipes(); 
    vector_t RRn = gerg_eos.Rgas_nodes(); 

    auto Z_pipes = gerg_eos.compute(temperature_pipes, pressure_pipes, x_pipes);
    auto Z_nodes = gerg_eos.compute(temperature_nodes, pressure_nodes, x_nodes);

    vector_t c2_pipes = Z_pipes.array() * RRp.array() * temperature_pipes.array(); 
    vector_t c2_nodes = Z_nodes.array() * RRn.array() * temperature_nodes.array(); 

    std::cout << __FILE__ << std::endl; 
    bool c2n_pass = verify_test("Speed of sound at nodes", c2_nodes, ref_c2_nodes);
    bool c2p_pass = verify_test("Speed of sound in pipes", c2_pipes, ref_c2_pipes);
          

    return !(c2n_pass & c2p_pass);
}


int main()
{
    bool pass_RR = RR_matlab();
    bool pass_matlab = gerg_matlab();
    bool pass_aga8cd = gerg_aga8code();

    
    return !(pass_matlab && pass_RR && pass_aga8cd);
}