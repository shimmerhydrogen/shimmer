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
#include "../src/viscosity.h"
using triple_t = std::array<double, 3>;

using namespace shimmer;

size_t num_steps = 25;
size_t num_inlet = 1;
size_t num_pipes = 15;
size_t num_nodes = 13;


enum station_type
{
    INLET,
    OUTLET,
    JUNCTION,
};


std::vector<station_type>
make_stations_type_vector()
{
    vector_t outlet_nodes(9);
    outlet_nodes << 1, 2, 5, 6, 8, 9, 10, 11, 12;

    vector_t inlet_nodes(1);
    inlet_nodes << 0; 

    vector_t junction_nodes(3);
    junction_nodes << 3, 4, 7; 

    std::vector<station_type> vec(num_nodes); 
    for(size_t i = 0; i < outlet_nodes.size(); i++)
        vec[outlet_nodes[i]] = station_type::OUTLET;            
    for(size_t i = 0; i < inlet_nodes.size(); i++)
        vec[inlet_nodes[i]] = station_type::INLET;            
    for(size_t i = 0; i < junction_nodes.size(); i++)
        vec[junction_nodes[i]] = station_type::JUNCTION;            
    
    return vec;
}


std::pair<vector_t, matrix_t>
read_boundary_data()
{
    vector_t Pset(num_steps);
    Pset << 70.000000000000000,  69.416666666666671,  68.833333333333343,
            68.250000000000000,  67.666666666666671,  67.083333333333343,
            66.500000000000000,  67.083333333333329,  67.666666666666657,
            68.249999999999986,  68.833333333333329,  69.416666666666657,
            70.000000000000000,  72.333333333333329,  74.666666666666657,
            76.999999999999972,  79.333333333333300,  81.666666666666629,
            84.000000000000000,  81.666666666666657,  79.333333333333329,
            77.000000000000000,  74.666666666666686,  72.333333333333357,
            70.000000000000000;
    Pset *=1E5;
    //---------------------------------------------------------------
    matrix_t Gsnam(num_steps, num_nodes);
    std::ifstream ifs("../unit_tests/gsnam.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam.txt not open" << std::endl;
        exit(1);
    }

    for (size_t icol = 0; icol < num_nodes; icol++)
        for (size_t irow = 0; irow < num_steps; irow++)
            ifs >>  Gsnam(irow, icol);
    ifs.close();

    return std::make_pair(Pset, Gsnam);
}

void
make_init_graph(infrastructure_graph& g)
{
    //---------------------------------------------------------------
    auto [Pset, Gsnam] = read_boundary_data();
    //---------------------------------------------------------------
    auto station_type_vec = make_stations_type_vector();
    //---------------------------------------------------------------
    std::vector<std::unique_ptr<station>> stations(num_nodes);

    for(size_t i = 0 ; i < num_nodes; i++)
    {   
        switch(station_type_vec[i])
        {
            case(station_type::INLET):
            {
                auto sin = make_inlet(Pset);
                stations[i] = std::make_unique<inlet_station>();
                stations[i]->set_state(sin);
                break;
            }
            case(station_type::OUTLET):
            {
                auto sout = make_outlet(Gsnam.col(i));
                stations[i] = std::make_unique<outlet_station>();
                stations[i]->set_state(sout);
                break;
            }
            case(station_type::JUNCTION):
                stations[i] = std::make_unique<junction>();
                break;
            default:
                throw std::invalid_argument("Station type not found");    
        }
    }
    //---------------------------------------------------------------

    std::vector<vertex_descriptor> vds;

    auto add_vertex = [&](vertex_properties&& vp) 
    {
        auto v = boost::add_vertex(g);
        g[v] = std::move(vp);
        return v;
    };

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( add_vertex({"station 0",  0, 70.000000000,-230 ,0.,x , std::move(stations[0])}));
    vds.push_back( add_vertex({"station 1",  1, 70.000000000,  20 ,0.,x , std::move(stations[1])}) );
    vds.push_back( add_vertex({"station 2",  2, 69.300000000,  25 ,0.,x , std::move(stations[2])}) );
    vds.push_back( add_vertex({"station 3",  3, 69.300000000,   0 ,0.,x , std::move(stations[3])}) );
    vds.push_back( add_vertex({"station 4",  4, 68.607000000,   0 ,0.,x , std::move(stations[4])}) );
    vds.push_back( add_vertex({"station 5",  5, 67.920930000,  20 ,0.,x , std::move(stations[5])}) );
    vds.push_back( add_vertex({"station 6",  6, 67.241720700,  30 ,0.,x , std::move(stations[6])}) );
    vds.push_back( add_vertex({"station 7",  7, 67.920930000,   0 ,0.,x , std::move(stations[7])}) );
    vds.push_back( add_vertex({"station 8",  8, 67.241720700,  50 ,0.,x , std::move(stations[8])}) );
    vds.push_back( add_vertex({"station 9",  9, 67.241720700,  20 ,0.,x , std::move(stations[9])}) );
    vds.push_back( add_vertex({"station 10",10, 66.569303493,  15 ,0.,x , std::move(stations[10])}) );
    vds.push_back( add_vertex({"station 11",11, 66.569303493,  40 ,0.,x , std::move(stations[11])}) );
    vds.push_back( add_vertex({"station 12",12, 67.241720700,  10 ,0.,x , std::move(stations[12])}) );

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



variable
make_guess_steady(size_t num_nodes, size_t num_pipes)
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);
    Gguess.setConstant(50);
    Pguess <<   70.000000000000000, 70.000000000000000, 69.299999999999997,
                69.299999999999997, 68.606999999999999, 67.920929999999998,
                67.241720700000002, 67.920929999999998, 67.241720700000002,
                67.241720700000002, 66.569303493000007, 66.569303493000007,
                67.241720700000002;

    //---------------------------------------------------------------
    // Read L from Gsnam / equal than in unsteady
    vector_t Lguess(num_nodes);
    std::ifstream ifs("../unit_tests/gsnam.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam.txt not open" << std::endl;
        exit(1);
    }

    for (size_t icol = 0; icol < num_nodes; icol++)
        ifs >>  Lguess(icol);
    ifs.close();


    return variable(Pguess, Gguess, Lguess); 
}



variable
make_guess_unsteady(size_t num_nodes, size_t num_pipes)
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);
    Pguess << 70.000000000000000, 66.142581067404251, 66.329787615521283,
              67.765953321847363, 63.706915015080810, 63.514270630072645,
              62.091515901977836, 61.486931580426337, 60.701860284128337,
              63.622896665695691, 55.649766619875237, 57.638531011806919,
              41.742844776429060;
    /*Gguess <<  239.0,    -25.0,  -94.0078045237288,   44.0078045237288,
              69.7886210532947,   20.0, 1.7376816440997, -92.0587439329238,
             -32.0587439329238,  -30.0, -10.0,-62.5, -38.7035744229765,
              16.5, 75.2035744229765;
    */ 
    Gguess << 2.390000000000000,  -0.250000000000000,  -0.940078045237288,
              0.440078045237288,   0.697886210532947,   0.200000000000000,
              0.017376816440997,  -0.920587439329238,  -0.320587439329238,
             -0.300000000000000,  -0.100000000000000,  -0.625000000000000,
             -0.387035744229765,   0.165000000000000,   0.752035744229765;    
    Pguess *= 1e5;
    Gguess *=1e2;
    //---------------------------------------------------------------
    // Read L from Gsnam
    vector_t Lguess(num_nodes);
    Lguess << -239, 25, 25, 0, 0,  20, 30, 0, 62.5, 20, 16.5,30, 10;

    #if 0
    std::ifstream ifs("../unit_tests/gsnam.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam.txt not open" << std::endl;
        exit(1);
    }

    for (size_t icol = 0; icol < num_nodes; icol++)
        ifs >>  Lguess(icol);
    ifs.close();
    #endif
    return variable(Pguess, Gguess, Lguess); 
}


std::pair<vector_t, matrix_t>
make_bnd_cond(size_t num_nodes, size_t num_pipes, size_t num_steps)
{
    size_t num_outlet = 9;    

    vector_t outlet_nodes(num_outlet);
    outlet_nodes << 1, 2, 5, 6, 8, 9, 10, 11, 12;

    vector_t Pset(num_steps);
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
    //---------------------------------------------------------------
    matrix_t Gsnam(num_steps, num_nodes);
    std::ifstream ifs("../unit_tests/gsnam.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam.txt not open" << std::endl;
        exit(1);
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

    return std::make_pair(Pset, flux_ext);  
}


std::pair<std::vector<double>, std::vector<double>>
make_reference(size_t num_nodes, size_t num_pipes)
{
    std::vector<double> ref_sol_unsteady = {7000000, 
        6934772.559942949, 6951327.384841953, 
        6994631.435108416, 6862997.00400152, 
        6846719.528087921, 6763531.784680299, 
        6756996.296753552, 6718661.845875668, 
        6898350.420853097, 6215288.565656994, 
        5636790.522900792, 3959225.464465233, 
        37.24788747552144, -24.06761888895476, 
       -52.4176615374912,   25.9037384225945, 
        40.34447212121921,  19.10036619989154, 
         1.223796919133901,-75.36356796841562, 
       -10.34670526035165, -49.06902068560661, 
       -12.11458220852576, -45.36865478624177,
       -32.85117746725343,  16.09781416668034, 
        36.45784656739564, -15.11186799782956, 
        25, 20,0, 0,20,30,0, 50,20, 16.5,50,12.5};


    variable var = make_guess_unsteady(num_nodes, num_pipes);
    vector_t vec = var.make_vector();

    std::vector<double> ref_sol_steady(num_nodes*2 + num_pipes);
    for(size_t i = 0; i< vec.size(); i++)
        ref_sol_steady[i] = vec(i);

    return std::make_pair(ref_sol_steady, ref_sol_unsteady);
}

int main()
{
    size_t num_steps = 25;
    size_t num_inlet = 1;
    size_t num_pipes = 15;
    size_t num_nodes = 13;

    size_t num_bcnd = num_nodes;
    size_t system_size = num_nodes + num_pipes;

    double tf = 25*3600;
    double dt = 3600;
    double temperature = 293.15;
    double tol = 1e-4;

    double tol_std = 1e-14; 
    double dt_std = 1;
    size_t MAX_ITER=1500;
  
    vector_t inlet_nodes(num_inlet);
    inlet_nodes << 0; 

    variable guess_unstd = make_guess_unsteady(num_nodes, num_pipes);
    variable guess_std   = make_guess_steady(num_nodes, num_pipes);
    auto [Pset, flux_ext] = make_bnd_cond(num_nodes, num_pipes, num_steps);
    //---------------------------------------------------------------

    infrastructure_graph graph;
    make_init_graph(graph);

    auto [y_nodes, y_pipes] = make_mass_fraction(num_nodes, graph);

    using time_solver_t = time_solver<papay, viscosity_type::Constant>; 

    
    //time_solver_t ts0(graph, temperature, Pset, flux_ext, inlet_nodes);
    //ts0.initialization(guess_std, dt_std, tol_std, y_nodes, y_pipes);    
    //auto sol_std =  ts0.guess();
    


    time_solver_t ts1(graph, temperature, Pset, flux_ext, inlet_nodes);
    ts1.set_initialization(guess_unstd);    
    ts1.advance(dt, num_steps, tol, y_nodes, y_pipes);
    auto sol_set_unstd  = ts1.solution();

    
    time_solver_t ts2(graph, temperature, Pset, flux_ext, inlet_nodes);
    ts2.initialization(guess_std, dt_std, tol_std, y_nodes, y_pipes);    
    ts2.advance(dt, num_steps, tol, y_nodes, y_pipes);
    auto sol_init_unstd = ts2.solution();
    
    //---------------------------------------------------------------
    std::cout << __FILE__ << std::endl; 
    auto [ref_std, ref_unstd] = make_reference(num_nodes, num_pipes);
    //bool pass0 = verify_test("time initialization", sol_std,  ref_std); 
    bool pass1 = verify_test("time solver with given init", sol_set_unstd, ref_unstd); 
    bool pass2 = verify_test("time solver computing  init", sol_init_unstd, ref_unstd); 

    return !(pass1 && pass2);
}
