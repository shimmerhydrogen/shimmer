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

#include "sol/sol.hpp"

#include "solver/incidence_matrix.h"
#include "solver/conservation_matrices.h"
#include "solver/fluid_solver.h"
#include "solver/time_solver.h"
#include "infra/infrastructure.h"
#include "errors.h"
struct config {
    std::string     database;
    size_t          steps;
    double          dt_std;
    double          dt;
    double          temperature;
    double          tol_std;
    double          tol;

    config();
};

config::config()
{
    database = "";
    steps = 1;
    dt_std = 1;
    dt = 3600;
    temperature = 293.15;
    tol_std = 1e-14;
    tol = 1e-4;
}

static void
register_usertypes(sol::state& lua)
{
    sol::usertype<config> config_ut = lua.new_usertype<config>(
        "config"
    );
    
    config_ut["database"] = &config::database;
    config_ut["steps"] = &config::steps;
    config_ut["dt_std"] = &config::dt_std;
    config_ut["dt"] = &config::dt;
    config_ut["temperature"] = &config::temperature;
    config_ut["tol_std"] = &config::tol_std;
    config_ut["tol"] = &config::tol;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Please specify configuration file name" << std::endl;
        return 1;
    }

    sol::state lua;
    lua.open_libraries(
        sol::lib::base, sol::lib::math, sol::lib::io,
        sol::lib::table, sol::lib::string
        );

    register_usertypes(lua);    

    config cfg;
    lua["config"] = &cfg;

    try {
        lua.safe_script_file(argv[1]);
    }
    catch (const sol::error& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    shimmer::infrastructure infra;

    int err = shimmer::load(cfg.database, infra);
    if (err != SHIMMER_SUCCESS) {
        std::cout << "Problem detected while loading DB" << std::endl;
        return 1;
    }

    shimmer::variable guess = initial_guess(infra);

    /* BEGIN GAS MASS FRACTIONS */
    int nstations = num_stations(infra);
    shimmer::matrix_t y_nodes = shimmer::matrix_t::Zero(nstations, NUM_GASES);
    for (size_t i = 0; i < infra.mass_fractions.size(); i++) {
        const auto& mf = infra.mass_fractions[i];
        assert(mf.i_snum < nstations);
        shimmer::vector_t y = shimmer::vector_t::Zero(NUM_GASES);
        std::copy(mf.fractions.begin(), mf.fractions.end(), y.begin());
        y_nodes.row(i) = y;
    }

    shimmer::incidence inc(infra.graph);
    shimmer::matrix_t y_pipes = inc.matrix_in().transpose() * y_nodes;   
    /* END GAS MASS FRACTIONS */

    using time_solver_t = shimmer::time_solver<shimmer::papay,
        shimmer::viscosity_type::Constant>;

    time_solver_t ts1(infra.graph, cfg.temperature);
    ts1.initialization(guess, cfg.dt_std, cfg.tol_std, y_nodes, y_pipes);  
    ts1.advance(cfg.dt, cfg.steps, cfg.tol, y_nodes, y_pipes);
    auto sol  = ts1.solution();
    std::cout << sol << std::endl;

    return 0;
}