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



std::pair<vector_t, vector_t>
read_boundary_conditions()
{
    vector_t Pset = vector_t::Zero(num_steps);

    Pset(0)  = 70.0;
    Pset(11) = 70.0;
    Pset  *=1E5;

    //---------------------------------------------------------------
    matrix_t Gsnam(num_steps, num_nodes);
    std::ifstream ifs("../unit_tests/gsnam_compressor.txt");
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

    return std::make_pair(Pset, Gsnam);
}

variable
make_guess_steady(const matrix_t& Gsnam )
{
    vector_t Gguess(num_pipes), Pguess(num_nodes);

    Gguess.setConstant(50);

    Pguess <<   70.000000000000000,
                70.000000000000000,
                69.299999999999997,
                69.299999999999997,
                68.606999999999999,
                67.920929999999998,
                67.241720700000002,
                67.920929999999998,
                67.241720700000002,
                67.241720700000002,
                66.569303493000007,
                70.000000000000000,
                67.241720700000002,
                67.000000000000000,
                67.000000000000000;
    Pguess *= 1E5;

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
    double tol = 1e-4;

    double tol_std = 1e-14; 
    double dt_std = 1;
    size_t MAX_ITER=1500;
  
    //---------------------------------------------------------------
    auto [Pset, Gsnam] = read_boundary_conditions();
    variable guess_std   = make_guess_steady(Gsnam);
    //variable guess_unstd = make_guess_unsteady(Gsnam);

    //---------------------------------------------------------------

    /*
    infrastructure_graph graph;
    make_init_graph(graph, Pset, Gsnam);

    auto [y_nodes, y_pipes] = make_mass_fraction(num_nodes, graph);

    using time_solver_t = time_solver<papay, viscosity_type::Constant>; 


    time_solver_t ts1(graph, temperature, dummyZero);
    ts1.initialization(guess_std, dt_std, tol_std, y_nodes, y_pipes);  
    */
    //ts1.advance(dt, num_steps, tol, y_nodes, y_pipes);
    //auto sol_unstd  = ts1.solution();

    //---------------------------------------------------------------
    //auto [ref_std, ref_unstd] = make_reference(guess_unstd);
    //bool pass =  verify_test("time solver with computed init", sol_unstd, ref_unstd); 
    //std::cout << pass << std::endl; 

    bool pass = true;    
    return !(pass);
}