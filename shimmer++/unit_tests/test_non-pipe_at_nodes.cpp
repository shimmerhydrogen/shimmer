/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it  
 */

#include <immintrin.h>
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
size_t num_pipes = 15;
size_t num_nodes = 13;
size_t num_outlet = 5;    

enum station_type
{
    REMI_WO_BACKFLOW,
    INJ_W_PRESS_CONTROL,
    OUTLET,
    JUNCTION,
    CONSUMPTION_WO_PRESS,
};


std::vector<station_type>
make_stations_type_vector()
{
    return std::vector<station_type>{
        station_type::REMI_WO_BACKFLOW,
        station_type::OUTLET,
        station_type::JUNCTION,
        station_type::JUNCTION,
        station_type::JUNCTION,
        station_type::CONSUMPTION_WO_PRESS,
        station_type::JUNCTION,
        station_type::JUNCTION,
        station_type::CONSUMPTION_WO_PRESS,
        station_type::JUNCTION,
        station_type::CONSUMPTION_WO_PRESS,
        station_type::INJ_W_PRESS_CONTROL,
        station_type::OUTLET
    };
}


vector_t
read_press_boundary_conditions()
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
    Pset  *=1E5;
    return Pset;
}


std::pair<matrix_t, matrix_t>
read_flux_boundary_conditions()
{
    vector_t outlet_nodes(num_outlet);
    outlet_nodes << 1,5,8,10,12;

    //---------------------------------------------------------------
    matrix_t Gsnam(num_steps, num_nodes);
    std::ifstream ifs("../unit_tests/gsnam_non_pipe.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam_non_pipe.txt not open" << std::endl;
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

    return std::make_pair(flux_ext, Gsnam);  
}



void
make_init_graph(infrastructure_graph& g, const vector_t& Pset, const matrix_t& Gsnam)
{
    //---------------------------------------------------------------
    double factor = 1.0;//0.85;
    //---------------------------------------------------------------
    std::vector<pair_input_t> user_constraints = {
                        std::make_pair(P_GREATER_EQUAL, 60E5),
                        std::make_pair(P_LOWER_EQUAL,   80E5),
                        std::make_pair(L_GREATER_EQUAL, -300),
                        std::make_pair(L_LOWER_EQUAL,    -10)};
    //---------------------------------------------------------------
    auto station_type_vec = make_stations_type_vector();
    //---------------------------------------------------------------
    std::vector<std::unique_ptr<station>> stations(num_nodes);

    for(size_t i = 0 ; i < num_nodes; i++)
    {   
        assert(i < station_type_vec.size());
        assert(i < stations.size());

        switch(station_type_vec[i])
        {
            case(station_type::JUNCTION):
                stations[i] = std::make_unique<junction>();
                break;
            case(station_type::REMI_WO_BACKFLOW):
            {
                auto remi = make_remi_wo_backflow(Pset, user_constraints,
                                                     user_constraints);
                stations[i] = std::make_unique<multiple_states_station>(remi);                
                break;
            }
            case(station_type::INJ_W_PRESS_CONTROL):
            {
                std::cout<< " INJ_W_PRESS: Lset"<< std::endl;
                for(size_t j = 0; j < Gsnam.rows(); j++)
                    std::cout << Gsnam(j,i) << std::endl;     

                auto inj_station = make_inj_w_pressure(factor, 7500000.0, Gsnam.col(i),
                                              user_constraints,
                                              user_constraints);
                stations[i] = std::make_unique<multiple_states_station>(inj_station);
                break;
            }
            case(station_type::CONSUMPTION_WO_PRESS):
            {
                auto consumption = make_consumption_wo_press(Gsnam.col(i), user_constraints);
                stations[i] = std::make_unique<one_state_station>(consumption);
                break;

            }
            case(station_type::OUTLET):
            {
                auto exit_station = make_outlet(Gsnam.col(i));
                stations[i] = std::make_unique<one_state_station>(exit_station);
                break;
            }
            default:
                throw std::invalid_argument("Station type not found");    
        }
    }
    //---------------------------------------------------------------

    std::vector<vertex_descriptor> vds;

    auto add_vertex = [&](vertex_properties&& vp, const vector_t& x_in, size_t i ) 
    {
        vp.gas_mixture = x_in;
        vp.node_station = std::move(stations[i]);

        auto v = boost::add_vertex(g);
        g[v] = std::move(vp);

        return v;
    };

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( add_vertex( vertex_properties("station 0",  0, 70.000000000, -75 ,0.),x , 0));
    vds.push_back( add_vertex( vertex_properties("station 1",  1, 70.000000000,  20 ,0.),x , 1) );
    vds.push_back( add_vertex( vertex_properties("station 2",  2, 69.300000000,   0 ,0.),x , 2) );
    vds.push_back( add_vertex( vertex_properties("station 3",  3, 69.300000000,   0 ,0.),x , 3) );
    vds.push_back( add_vertex( vertex_properties("station 4",  4, 68.607000000,   0 ,0.),x , 4) );
    vds.push_back( add_vertex( vertex_properties("station 5",  5, 67.920930000,  20 ,0.),x , 5) );
    vds.push_back( add_vertex( vertex_properties("station 6",  6, 67.241720700,   0 ,0.),x , 6) );
    vds.push_back( add_vertex( vertex_properties("station 7",  7, 67.920930000,   0 ,0.),x , 7) );
    vds.push_back( add_vertex( vertex_properties("station 8",  8, 67.241720700,  50 ,0.),x , 8) );
    vds.push_back( add_vertex( vertex_properties("station 9",  9, 67.241720700,   0 ,0.),x , 9) );
    vds.push_back( add_vertex( vertex_properties("station 10",10, 66.569303493,  15 ,0.),x , 10) );
    vds.push_back( add_vertex( vertex_properties("station 11",11, 70.000000000, -40 ,0.),x , 11) );
    vds.push_back( add_vertex( vertex_properties("station 12",12, 67.241720700,  10 ,0.),x , 12) );

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
make_guess_steady()
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);
    Gguess.setConstant(50);
    Pguess <<   70.000000000000000,   70.000000000000000,
                69.299999999999997,   69.299999999999997,
                68.606999999999999,   67.920929999999998,
                67.241720700000002,   67.920929999999998,
                67.241720700000002,   67.241720700000002,
                66.569303493000007,   70.000000000000000,
                67.241720700000002;
    Pguess *= 1E5;

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
make_guess_unsteady(const matrix_t& Gsnam)
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);

    Pguess << 7.000000000000000,   6.952086053520747,   6.963715087477524,
              6.976886223836662,   6.953982373265999,   6.936503561559128,
              6.953818900684214,   6.862710129421968,   6.817393691883836,
              6.901609063012971,   6.314666242105152,   7.583706514628937,
              5.197436151747834;
    Pguess *= 1e6;

    Gguess <<  75.000000000000000, -20.000000000000000, -27.811646427934370,
                7.811646427934369,  15.684394159687432,  20.000000000000000,
                1.135918148126126,  -2.360122439495675, -42.360122439495676,
               40.000000000000000, -10.000000000000000, -50.000000000000000,
              -16.503959412378201,  15.000000000000000,  31.503959412378201;
    vector_t Lguess = Gsnam.row(1);

    return variable(Pguess, Gguess, Lguess); 
}


std::pair<std::vector<double>, std::vector<double>>
make_reference(variable& guess_unsteady)
{
    std::vector<double> ref_sol_unsteady = {6689825.904715838,
  6681553.365209631,
  6688466.308533575,
  6689698.724947638,
  6690387.967752758,
  6678818.878093707,
  6693392.727712549,
  6607500.072127393,
  6541462.432020029,
  6659051.37492536,
  6175980.946977876,
  7500000,
  5549672.589500213,
  4.339917747620353,
  -14.78723050483148,
  -7.528609364865357,
  -3.142642224954374,
  -2.179240479202206,
  15.68166036749584,
  1.052436381983676,
  11.37359735961809,
  -40.18692744914854,
  44.70630625247412,
  -8.036112103250305,
  -59.53594638078679,
  -18.78889725554646,
  13.40415996798448,
  18.88014903025844,
  0,
  15,
  0,
  0,
  0,
  16,
  0,
  0,
  62.5,
  0,
  13.5,
  -44.70630625247412,
  8
};


    vector_t vec = guess_unsteady.make_vector();

    std::vector<double> ref_sol_steady(num_nodes*2 + num_pipes);
    for(size_t i = 0; i< vec.size(); i++)
        ref_sol_steady[i] = vec(i);

    return std::make_pair(ref_sol_steady, ref_sol_unsteady);
}


int main()
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    size_t num_bcnd = num_nodes;
    size_t system_size = num_nodes + num_pipes;

    size_t num_steps = 7;
    double dt = 3600;
    double temperature = 293.15;
    double tol = 1e-4;

    double tol_std = 1e-14; 
    double dt_std = 1;
    size_t MAX_ITER=1500;
  
    vector_t Pset = read_press_boundary_conditions();
    auto [flux_ext, Gsnam] = read_flux_boundary_conditions();
    variable guess_unstd = make_guess_unsteady(Gsnam);
    variable guess_std   = make_guess_steady();

    //---------------------------------------------------------------

    infrastructure_graph graph;
    make_init_graph(graph, Pset, Gsnam);


    auto [y_nodes, y_pipes] = make_mass_fraction(num_nodes, graph);

    using time_solver_t = time_solver<papay, viscosity_type::Constant>; 

    #if 0
    {
    time_solver_t ts0(graph, temperature, flux_ext);
    ts0.set_initialization(guess_unstd);    
    ts0.advance(dt, num_steps, tol, y_nodes, y_pipes);
    auto sol_set_unstd  = ts0.solution();
    }
    #endif

    time_solver_t ts1(graph, temperature, flux_ext);
    ts1.initialization(guess_std, dt_std, tol_std, y_nodes, y_pipes);  
    ts1.advance(dt, num_steps, tol, y_nodes, y_pipes);
    auto sol_unstd  = ts1.solution();

    //---------------------------------------------------------------
    auto [ref_std, ref_unstd] = make_reference(guess_unstd);
    bool pass =  verify_test("time solver with computed init", sol_unstd, ref_unstd); 

    std::cout << pass << std::endl; 

    return !(pass);
}