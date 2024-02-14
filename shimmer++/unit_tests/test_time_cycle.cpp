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
#include "../src/time_solver.h"

using triple_t = std::array<double, 3>;

using namespace shimmer;


void
make_init_graph(infrastructure_graph& g)
{

    std::vector<vertex_descriptor> vds;

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( boost::add_vertex( { "station 0",  0, 70.000000000,-230, 0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 1",  1, 70.000000000,  20 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 2",  2, 69.300000000,  25 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 3",  3, 69.300000000,   0 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 4",  4, 68.607000000,   0 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 5",  5, 67.920930000,  20 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 6",  6, 67.241720700,  30 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 7",  7, 67.920930000,   0 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 8",  8, 67.241720700,  50 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 9",  9, 67.241720700,  20 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 10",10, 66.569303493,  15 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 11",11, 66.569303493,  40 ,0.,x }, g) );
    vds.push_back( boost::add_vertex( { "station 12",12, 67.241720700,  10 ,0.,x }, g) );


    using eprop_t = edge_properties;

    edge_properties ep0  = {edge_type::pipe, 0,  80000,	1.2,	1.20E-05};
    edge_properties ep1  = {edge_type::pipe, 1,  16000,	0.6,	1.20E-05};
    edge_properties ep2  = {edge_type::pipe, 2,  40000,	0.8,	1.20E-05};
    edge_properties ep3  = {edge_type::pipe, 3, 160000,	0.7,	1.20E-05};
    edge_properties ep4  = {edge_type::pipe, 4, 200000,	0.8,	1.20E-05};
    edge_properties ep5  = {edge_type::pipe, 5,  24000,	0.6,	1.20E-05};
    edge_properties ep6  = {edge_type::pipe, 6, 120000,	0.2,	1.20E-05};
    edge_properties ep7  = {edge_type::pipe, 7,  80000,	0.9,	1.20E-05};
    edge_properties ep8  = {edge_type::pipe, 8,  64000,	0.7,	1.20E-05};
    edge_properties ep9  = {edge_type::pipe, 9, 240000,	0.6,	1.20E-05};
    edge_properties ep10 = {edge_type::pipe,10,  28000,	0.2,	1.20E-05};
    edge_properties ep11 = {edge_type::pipe,11,  80000,	0.9,	1.20E-05};
    edge_properties ep12 = {edge_type::pipe,12, 160000,	0.7,	1.20E-05};
    edge_properties ep13 = {edge_type::pipe,13,  40000,	0.3,	1.20E-05};
    edge_properties ep14 = {edge_type::pipe,14, 320000,	0.9,	1.20E-05};


    boost::add_edge( vds[ 0], vds[ 3], ep0, g);
    boost::add_edge( vds[ 1], vds[ 2], ep1, g);
    boost::add_edge( vds[ 2], vds[ 3], ep2, g);
    boost::add_edge( vds[ 2], vds[ 4], ep3, g);
    boost::add_edge( vds[ 3], vds[ 4], ep4, g);
    boost::add_edge( vds[ 4], vds[ 5], ep5, g);
    boost::add_edge( vds[ 4], vds[ 7], ep6, g);
    boost::add_edge( vds[ 6], vds[ 4], ep7, g);
    boost::add_edge( vds[ 7], vds[ 6], ep8, g);
    boost::add_edge( vds[11], vds[ 6], ep9, g);
    boost::add_edge( vds[12], vds[ 7], ep10, g);
    boost::add_edge( vds[ 8], vds[ 7], ep11, g);
    boost::add_edge( vds[ 7], vds[ 9], ep12, g);
    boost::add_edge( vds[ 9], vds[10], ep13, g);
    boost::add_edge( vds[ 3], vds[ 9], ep14, g);

}


std::pair<matrix_t, matrix_t> 
make_mass_fraction(size_t size, const infrastructure_graph& graph)
{
    incidence inc(graph);

    matrix_t y_nodes(size, 21);
    y_nodes.col(0).setConstant(1);

    matrix_t y_pipes = inc.matrix_in().transpose() * y_nodes;    


    return  std::make_pair(y_nodes, y_pipes); 
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
    size_t num_steps = 25;

    size_t num_inlet = 1;
    size_t num_outlet = 9;    

    size_t num_pipes = 15;
    size_t num_nodes = 13;

    size_t num_bcnd = num_nodes;
    size_t system_size = num_nodes + num_pipes;

    double tf = 25*3600;
    double dt = 3600;
    double temperature = 293.15;
    double tolerance = 1e-6;

    vector_t Pguess(num_nodes), Gguess(num_pipes), Lguess(num_nodes);        
    vector_t Pset(num_steps), inlet_nodes(num_inlet), outlet_nodes(num_outlet);

    inlet_nodes << 0; 
    outlet_nodes << 1, 2, 5, 6, 8, 9, 10, 11, 12;
    Pguess << 70.000000000000000, 66.142581067404251, 66.329787615521283,
              67.765953321847363, 63.706915015080810, 63.514270630072645,
              62.091515901977836, 61.486931580426337, 60.701860284128337,
              63.622896665695691, 55.649766619875237,  57.638531011806919,
              41.742844776429060;
    Gguess <<  239.0, -25.0, -94.0078045237288, 44.0078045237288,
              69.7886210532947,   20.0, 1.7376816440997, -92.0587439329238,
             -32.0587439329238, -30.0, -10.0, -62.5, -38.7035744229765,
              16.5, 75.2035744229765;
    Pset << 70.000000000000000,  69.416666666666671,  68.833333333333343,
            68.250000000000000,  67.666666666666671,  67.083333333333343,
            66.500000000000000,  67.083333333333329,  67.666666666666657,
            68.249999999999986,  68.833333333333329,  69.416666666666657,
            70.000000000000000,  72.333333333333329,  74.666666666666657,
            76.999999999999972,  79.333333333333300,  81.666666666666629,
            84.000000000000000,  81.666666666666657,  79.333333333333329,
            77.000000000000000,  74.666666666666686,  72.333333333333357,
            70.000000000000000;
    Pset  *=1E5;
    Pguess *= 1e5;

    matrix_t Gsnam(num_steps, num_nodes);
    std::ifstream ifs("../unit_tests/gsnam.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam.txt not open" << std::endl;
        return 1;
    }

    for (size_t icol = 0; icol < num_nodes; icol++)
        for (size_t irow = 0; irow < num_steps; irow++)
            ifs >>  Gsnam(irow, icol);
    ifs.close();


    //---------------------------------------------------------------
    matrix_t flux_ext = matrix_t::Zero(num_steps, num_nodes);
    for(size_t i= 0 ; i < outlet_nodes.size(); i++)
    {   size_t idx = outlet_nodes(i);
        flux_ext.col(idx) = Gsnam.col(idx);     
    }

    Lguess =  Gsnam.row(0);
    //---------------------------------------------------------------

    infrastructure_graph graph;
    make_init_graph(graph);

    auto [y_nodes, y_pipes] = make_mass_fraction(num_nodes, graph);

    time_solver<papay> ts(graph, temperature, Pset, flux_ext, inlet_nodes);
    ts.initialization(Pguess, Gguess, Lguess);
    ts.advance(dt, num_steps, tolerance, graph, y_nodes, y_pipes);

    vector_t sol(num_bcnd + num_pipes + num_nodes);
    //sol.head(num_nodes) = var.pressure;
    //sol.segment(num_nodes, num_pipes) = var.flux;
    //sol.tail(num_bcnd) = var.L_rate;

    bool pass = false;//verify_test("Test fluid-dynamic solver", sol, ref_sol); 


    return !pass;
}
