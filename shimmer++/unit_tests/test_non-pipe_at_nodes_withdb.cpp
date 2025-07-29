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

#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "solver/conservation_matrices.h"
#include "solver/fluid_solver.h"
#include "solver/time_solver.h"
#include "solver/viscosity.h"
#include "verify_test.h"
#include "infra/infrastructure.h"

#include "sqlite/sqlite.hpp"
#include "errors.h"

using triple_t = std::array<double, 3>;

using namespace shimmer;

size_t num_steps = 25;
size_t num_pipes = 15;
size_t num_nodes = 13;
size_t num_outlet = 5;    

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



variable
make_guess_steady()
{
    vector_t Gguess(::num_pipes), Pguess(num_nodes);
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
    vector_t Gguess(::num_pipes), Pguess(num_nodes);

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
    vector_t Lguess = Gsnam.row(0);
    std::cout << "LGuess: " << Lguess.transpose() << std::endl;

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

    std::vector<double> ref_sol_steady(num_nodes*2 + ::num_pipes);
    for(size_t i = 0; i< vec.size(); i++)
        ref_sol_steady[i] = vec(i);

    return std::make_pair(ref_sol_steady, ref_sol_unsteady);
}





int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Please specify database file name" << std::endl;
        return 1;
    }

    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    

    size_t num_steps = 7;
    double dt = 3600;
    double temperature = 293.15;
    double tol = 1e-4;

    double tol_std = 1e-14; 
    double dt_std = 1;
    size_t MAX_ITER=1500;
  
    auto [flux_ext, Gsnam] = read_flux_boundary_conditions();
    variable guess_unstd = make_guess_unsteady(Gsnam);
    //variable guess_std   = make_guess_steady();
    //---------------------------------------------------------------

    shimmer::infrastructure infra;

    int err = shimmer::load(argv[1], infra);
    if (err != SHIMMER_SUCCESS) {
        std::cout << "Problem detected while loading DB" << std::endl;
        return 1;
    }

    size_t num_bcnd = num_stations(infra);
    size_t system_size = num_stations(infra) + shimmer::num_pipes(infra);
    variable guess_std = initial_guess(infra);

    using time_solver_t = time_solver<papay, viscosity_type::Constant>; 

    #if 0
    {
    time_solver_t ts0(graph, temperature, flux_ext);
    ts0.set_initialization(guess_unstd);    
    ts0.advance(dt, num_steps, tol);
    auto sol_set_unstd  = ts0.solution();
    }
    #endif

    time_solver_t ts1(infra.graph, temperature);
    ts1.initialization(guess_std, dt_std, tol_std);  
    ts1.advance(dt, num_steps, tol);
    auto sol_unstd  = ts1.solution();
    auto sol_std  = ts1.guess();
    //---------------------------------------------------------------
    auto [ref_std, ref_unstd] = make_reference(guess_unstd);
    auto xx = guess_unstd.make_vector();
    bool pass =  verify_test("time solver with computed init", sol_std, {xx.begin(), xx.end()} ); 

    std::cout << pass << std::endl; 

    return !(pass);
}