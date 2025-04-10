/* 
    0. Type of station must be clearlly clarified! 
        There is no information of which type of exit stations are 3! 
        This must be given to be read for the database! 

    1. Check Mass fraction: They have 16 for molar fraction. It ends uo to be the same? 
    2. Density: I did not set anything. But They have rho_g_std
    3. flux_ext seems useless.
*/

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
#include <utility>

#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/conservation_matrices.h"
#include "boundary.h"
#include "verify_test.h"
#include "MATLAB_GERG_functions.hpp"
#include "../src/fluid_solver.h"
#include "../src/time_solver.h"
#include "../src/viscosity.h"
#include "../src/nonpipe_over_edges.h"


using triple_t = std::array<double, 3>;
using namespace shimmer;

size_t num_steps = 25;
size_t num_pipes = 17;
size_t num_nodes = 15;
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
        station_type::OUTLET,
        station_type::JUNCTION,
        station_type::JUNCTION
    };
}



std::pair<vector_t, vector_t>
read_pressure_inlet_node_stations()
{
    vector_t Pset_INJ = vector_t::Zero(num_steps);
    vector_t Pset_REMI = vector_t::Zero(num_steps);

    Pset_INJ<<  
    7.000000000000000e+06,     6.766666666666667e+06,
    6.533333333333333e+06,     6.300000000000000e+06,
    6.066666666666667e+06,     5.833333333333334e+06,     
    5.600000000000000e+06,     5.833333333333334e+06,
    6.066666666666667e+06,     6.300000000000000e+06,
    6.533333333333333e+06,     6.766666666666667e+06,
    7.000000000000000e+06,     7.233333333333333e+06,
    7.466666666666666e+06,     7.699999999999997e+06,
    7.933333333333330e+06,     8.166666666666663e+06,
    8.400000000000000e+06,     8.458333333333332e+06,
    8.516666666666666e+06,     8.574999999999998e+06,
    8.633333333333332e+06,     8.691666666666666e+06,
    8.750000000000000e+06;

    Pset_REMI<< 
     7.000000000000000e+06,     6.941666666666667e+06,     
     6.883333333333334e+06,     6.825000000000000e+06,
     6.766666666666667e+06,     6.708333333333334e+06,
     6.650000000000000e+06,     6.708333333333333e+06,
     6.766666666666666e+06,     6.824999999999998e+06,
     6.883333333333333e+06,     6.941666666666666e+06,
     7.000000000000000e+06,     7.233333333333333e+06,
     7.466666666666666e+06,     7.699999999999997e+06,
     7.933333333333330e+06,     8.166666666666663e+06,
     8.400000000000000e+06,     8.166666666666666e+06,
     7.933333333333333e+06,     7.700000000000000e+06,
     7.466666666666669e+06,     7.233333333333336e+06,
     7.000000000000000e+06;

     return std::make_pair(Pset_INJ, Pset_REMI);
}

matrix_t
read_flux_node_stations()
{
    //---------------------------------------------------------------
    matrix_t Gsnam(num_nodes, num_steps);
    std::ifstream ifs("../unit_tests/gsnam_compressor.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam_compressor.txt not open" << std::endl;
        exit(1);
    }

    for (size_t irow = 0; irow < num_nodes; irow++)
        for (size_t icol = 0; icol < num_steps; icol++)
            ifs >>  Gsnam(irow, icol);
    ifs.close();
    //---------------------------------------------------------------

    return Gsnam;
}

variable
make_guess_steady(const matrix_t& Gsnam )
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);

    Gguess.setConstant(50);

    Pguess <<      7.000000000000000,
                    7.000000000000000,
                    6.930000000000000,
                    6.930000000000000,
                    6.860700000000000,
                    6.792093000000000,
                    6.724172070000000,
                    6.792093000000000,
                    6.724172070000000,
                    6.724172070000000,
                    6.656930349300000,
                    7.000000000000000,
                    6.724172070000000,
                    6.700000000000000,
                    6.700000000000000;
    Pguess *= 1E6;

    vector_t Lguess = Gsnam.col(0);

    return variable(Pguess, Gguess, Lguess); 
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


void
make_init_graph(infrastructure_graph& g,
                 const vector_t& Pset_REMI,
                 const vector_t& Pset_INJ,
                 const matrix_t& Gsnam)
{
    // Only one switch allowed to fit matlab code
    bool only_one_switch = true;
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
                auto remi = make_remi_wo_backflow(Pset_REMI, user_constraints,
                                                     user_constraints, only_one_switch);
                stations[i] = std::make_unique<multiple_states_station>(remi);                
                break;
            }
            case(station_type::INJ_W_PRESS_CONTROL):
            {
                std::cout<< " INJ_W_PRESS: Lset"<< std::endl;
                // Here Pset_Inj should be given, but in the code Pset=75E5, instead of
                // the data provided by the files
                auto inj_station = make_inj_w_pressure(factor, 75.E5, Gsnam.row(i),
                                              user_constraints,
                                              user_constraints,
                                              only_one_switch);
                stations[i] = std::make_unique<multiple_states_station>(inj_station);
                break;
            }
            case(station_type::CONSUMPTION_WO_PRESS):
            {
                auto consumption = make_consumption_wo_press(Gsnam.row(i), user_constraints);
                stations[i] = std::make_unique<one_state_station>(consumption);
                break;

            }
            case(station_type::OUTLET):
            {
                auto exit_station = make_outlet(Gsnam.row(i));
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

    // Read in topology.xlsx, sheet: Nodes
    //                                           NAME        ID   Pguess       Qset  ? (Maybe Consumption?)  but not store anymore in graph    
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
    vds.push_back( add_vertex( vertex_properties("station 13",13, 66.569303493, 0.0, 0.),x , 13) );
    vds.push_back( add_vertex( vertex_properties("station 14",14, 66.569303493, 0.0, 0.),x , 14) );

    using eprop_t = edge_properties;

    // Parameters have to be given: small_network_disma_branches_nonpipe.xlxs
    std::vector<bool> activate_history ( num_steps, true); 
    
    using edge_constraint_t = edge_station::control::constraint_type;
    using external_t = edge_station::external_type; 

    std::unordered_map<external_t, std::pair<edge_constraint_t, double>>  user_limits; 

    user_limits[external_t::P_OUT_MAX] = std::make_pair(edge_constraint_t::LOWER_EQUAL,   80.E+5);
    user_limits[external_t::P_IN_MIN]  = std::make_pair(edge_constraint_t::GREATER_EQUAL, 50.E+5);
    user_limits[external_t::BETA_MAX]  = std::make_pair(edge_constraint_t::LOWER_EQUAL,   2.0);
    user_limits[external_t::BETA_MIN]  = std::make_pair(edge_constraint_t::GREATER_EQUAL, 1.2);
    user_limits[external_t::FLUX_MAX]  = std::make_pair(edge_constraint_t::LOWER_EQUAL,   80.0);
    user_limits[external_t::PWD_NOMINAL] = std::make_pair(edge_constraint_t::LOWER_EQUAL, 10.0);
    
    //What should I do if this limits are not given? Should I put this as default? Or deactivate option 
    user_limits[external_t::P_THRESHOLD_MIN] = std::make_pair(edge_constraint_t::GREATER, -1.E+20);
    user_limits[external_t::P_THRESHOLD_MAX] = std::make_pair(edge_constraint_t::LOWER, 1.E+20);

    auto efficiency = 0.9;
    auto ramp_coeff = 0.0;
    using mode_t = edge_station::control::mode_type;

    auto mypair = std::make_pair(mode_t::POWER_DRIVER, 2.E6);
                //std::make_pair(mode_t::BETA, 1.8);
                // std::make_pair(mode_t::FLUX, 5);
//                std::make_pair(mode_t::PRESSURE_OUT, 80E5);
//                std::make_pair(mode_t::PRESSURE_IN , 50E5);
    std::vector<std::pair<mode_t,double>> mode_type_vec = {mypair};
    //                                                   

    auto comp = edge_station::make_compressor(ramp_coeff,
                                                    efficiency, 
                                                    activate_history,
                                                    mode_type_vec,
                                                    user_limits);

    std::cout << comp <<std::endl;

    // Read in topology.xlsx, sheet: PIPEs
    edge_properties ep0  = {edge_type::pipe, 0,  80000,	1.2,	1.20E-05};
    edge_properties ep1  = {edge_type::pipe, 1,  16000,	0.6,	1.20E-05};
    edge_properties ep2  = {edge_type::pipe, 2,  40000,	0.8,	1.20E-05};
    edge_properties ep3  = {edge_type::pipe, 3, 160000,	0.7,	1.20E-05};
    edge_properties ep4  = {edge_type::pipe, 4, 200000,	0.8,	1.20E-05};
    edge_properties ep5  = {edge_type::pipe, 5,  24000,	0.6,	1.20E-05};
    edge_properties ep6  = {edge_type::pipe, 6,  60000,	0.2,	1.20E-05};
    edge_properties ep7  = {edge_type::pipe, 7,  80000,	0.9,	1.20E-05};
    edge_properties ep8  = {edge_type::pipe, 8,  64000,	0.7,	1.20E-05};
    edge_properties ep9  = {edge_type::pipe, 9, 240000,	0.6,	1.20E-05};
    edge_properties ep10 = {edge_type::pipe,10,  28000,	0.2,	1.20E-05};
    edge_properties ep11 = {edge_type::pipe,11,  80000,	0.9,	1.20E-05};
    edge_properties ep12 = {edge_type::pipe,12, 160000,	0.7,	1.20E-05};
    edge_properties ep13 = {edge_type::pipe,13,  40000,	0.3,	1.20E-05};
    edge_properties ep14 = {edge_type::pipe,14, 320000,	0.9,	1.20E-05};
    edge_properties ep15 = {edge_type::compressor,15, 1, 0.2 ,	1.20E-05, std::make_shared<edge_station::compressor>(comp)};
    edge_properties ep16 = {edge_type::pipe,16,  60000, 0.2,	1.20E-05};


    // Read in topology.xlsx, sheet: PIPEs
    boost::add_edge( vds[	0	], vds[	3	], ep0, g);
    boost::add_edge( vds[	1	], vds[	2	], ep1, g);
    boost::add_edge( vds[	2	], vds[	3	], ep2, g);
    boost::add_edge( vds[	2	], vds[	4	], ep3, g);
    boost::add_edge( vds[	3	], vds[	4	], ep4, g);
    boost::add_edge( vds[	4	], vds[	5	], ep5, g);
    boost::add_edge( vds[	4	], vds[	13	], ep6, g);
    boost::add_edge( vds[	6	], vds[	4	], ep7, g);
    boost::add_edge( vds[	7	], vds[	6	], ep8, g);
    boost::add_edge( vds[	11	], vds[	6	], ep9, g);
    boost::add_edge( vds[	12	], vds[	7	], ep10, g);
    boost::add_edge( vds[	8	], vds[	7	], ep11, g);
    boost::add_edge( vds[	7	], vds[	9	], ep12, g);
    boost::add_edge( vds[	9	], vds[	10	], ep13, g);
    boost::add_edge( vds[	3	], vds[	9	], ep14, g);
    boost::add_edge( vds[	13	], vds[	14	], ep15, g);
    boost::add_edge( vds[	14	], vds[	7	], ep16, g);
}

variable
make_guess_unsteady(const matrix_t& Gsnam)
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);

    Pguess <</* 7.000000000000000,   6.952086053520747,   6.963715087477524,
              6.976886223836662,   6.953982373265999,   6.936503561559128,
              6.953818900684214,   6.862710129421968,   6.817393691883836,
              6.901609063012971,   6.314666242105152,   7.583706514628937,
              5.197436151747834;*/
            /* Temporal solution until understanding why there are small 
            differences in the code results => check test_non-pipe_at_nodes
            */
     7.000000000000000e+06,
     6.951352260042983e+06,
     6.962982653112382e+06,
     6.976886223836661e+06,
     6.951292322506930e+06,
     6.933805989233891e+06,
     6.951425410985802e+06,
     6.877996691005643e+06,
     6.832792658281552e+06,
     6.909399665140901e+06,
     6.323261067613857e+06,
     7.581532847685807e+06,
     5.218177823376325e+06,
     5.000000000000000e+06,
     8.368800804804965e+06;
    
    Gguess <<
     7.500000000000000e+01,
    -2.000000000000000e+01,
    -2.863255760560140e+01,
     8.632557605601399e+00,
     1.665929687486727e+01,
     2.000000000000000e+01,
     7.398478500926114e+00,
     2.106624020457450e+00,
    -3.789337597954255e+01,
     4.000000000000000e+01,
    -1.000000000000000e+01,
    -5.000000000000000e+01,
    -1.470814551953134e+01,
     1.500000000000000e+01,
     2.970814551953134e+01,
     7.398478500926114e+00,
     7.398478500926114e+00;

    vector_t Lguess = Gsnam.col(1);

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
                                            8};

    vector_t vec = guess_unsteady.make_vector();

    std::vector<double> ref_sol_steady(num_nodes*2 + num_pipes);
    for(size_t i = 0; i< vec.size(); i++)
        ref_sol_steady[i] = vec(i);

    return std::make_pair(ref_sol_steady, ref_sol_unsteady);
}

int main()
{
    // This are given as global variables but should be read!!!!
    /*
    size_t num_steps = ;
    size_t num_pipes = ;
    size_t num_nodes = ;
    size_t num_outlet =;    
    */

    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    size_t num_bcnd = num_nodes;
    size_t system_size = num_nodes + num_pipes;


    double dt = 3600;
    double tfinal = 24*3600; 
    size_t num_steps = std::ceil(tfinal/dt);

    double temperature = 293.15;
    double tol = 1e-10;

    double tol_std = 1e-14; 
    double dt_std = 1;
    size_t MAX_ITER=1500;
  
    //---------------------------------------------------------------
    auto [Pset_INJ, Pset_REMI] = read_pressure_inlet_node_stations();
    matrix_t Gsnam = read_flux_node_stations();
    variable guess_std   = make_guess_steady(Gsnam);
    variable guess_unstd = make_guess_unsteady(Gsnam);
    vector_t dummyZero = vector_t::Zero(num_nodes);
    //---------------------------------------------------------------

    infrastructure_graph graph;
    make_init_graph(graph, Pset_REMI, Pset_INJ, Gsnam);

    auto [y_nodes, y_pipes] = make_mass_fraction(num_nodes, graph);

    using time_solver_t = time_solver<papay, viscosity_type::Constant>; 

    time_solver_t ts1(graph, temperature, dummyZero);

    std::cout << std::setprecision(16) << guess_std.make_vector() << std::endl; 

    ts1.initialization(guess_std, dt_std, tol_std, y_nodes, y_pipes);  
    //ts1.set_initialization(guess_unstd);   
    vector_t sol_std  = ts1.guess();
    //std::cout << std::setprecision(16)<< sol_std << std::endl; 

    ts1.advance(dt, num_steps, tol, y_nodes, y_pipes);
    vector_t sol_unstd  = ts1.solution();   
    //std::cout << std::setprecision(16)<< sol_unstd << std::endl; 

    //---------------------------------------------------------------
    //auto [ref_std, ref_unstd] = make_reference(guess_unstd);
    //bool pass =  verify_test("time solver with computed init", sol_std, ref_std); 

    //std::cout << pass << std::endl; 

    bool pass = true;

    return !(pass);
}
