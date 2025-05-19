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

#define SOL_ALL_SAFETIES_ON 1

#include <iostream>
#include <unistd.h>

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

static int
parse_cmdline(int& argc, char**&argv, std::string& luafn, config& cfg)
{
    int opt;
    while ( (opt = getopt(argc, argv, "i:")) != -1) {
        switch (opt) {
        case 'i':
            cfg.database = optarg;
            break;
        default:
            std::cerr << "Error parsing cmdline" << std::endl;
            return -1;
            break;
        }
    }

    if (cfg.database == "") {
        std::cerr << "Database not specified (either with ";
        std::cerr << "-i flag or from config file)";
        std::cerr << std::endl;
        return -1;
    }

    return 0;
}

namespace shimmer {

int save_pressures(const std::string& db_filename, const infrastructure& infra,
    const matrix_t& pressures)
{
    assert(pressures.cols() == infra.s_i2u.size());

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_station_pressures VALUES (?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_station_pressures", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    
    for (int ts = 0; ts < pressures.rows(); ts++) {
        for (int i_snum = 0; i_snum < pressures.cols(); i_snum++) {
            rc = sqlite3_bind_int(stmt, 1, convert_i2u(infra.s_i2u, i_snum));
            rc = sqlite3_bind_int(stmt, 2, ts);
            rc = sqlite3_bind_double(stmt, 3, pressures(ts,i_snum));
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}

int save_flowrates(const std::string& db_filename, const infrastructure& infra,
    const matrix_t& flowrates)
{
    assert(flowrates.cols() == num_pipes(infra));

    auto edge_range = boost::edges(infra.graph);
    auto edge_begin = edge_range.first;
    auto edge_end = edge_range.second;

    sqlite3 *db = nullptr;
    int rc = sqlite3_open_v2(db_filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr);
    if(rc) {
        std::cerr << "Can't open database '" << db_filename << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO solution_pipe_flowrates VALUES (?, ?, ?, ?, ?)";

    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }
    
    rc = sqlite3_exec(db, "DELETE FROM solution_pipe_flowrates", nullptr, nullptr, nullptr);
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    
    for (int ts = 0; ts < flowrates.rows(); ts++) {
        for (int branch_num = 0; branch_num < flowrates.cols(); branch_num++) {
            auto edge_itor = std::next(edge_begin, branch_num);
            auto ep = infra.graph[*edge_itor];
            rc = sqlite3_bind_text(stmt, 1, ep.name.c_str(), -1, nullptr);
            rc = sqlite3_bind_int(stmt, 2, convert_i2u(infra.s_i2u, ep.i_sfrom));
            rc = sqlite3_bind_int(stmt, 3, convert_i2u(infra.s_i2u, ep.i_sto));
            rc = sqlite3_bind_int(stmt, 4, ts);
            rc = sqlite3_bind_double(stmt, 5, flowrates(ts, branch_num));
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
    }
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    sqlite3_close(db);
    return SHIMMER_SUCCESS;
}


} //namespace shimmer

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

    auto sresult = lua.safe_script_file(argv[1], sol::script_pass_on_error);
    if (!sresult.valid()) {
        sol::error err = sresult;
        std::cerr << "Solver configuration: " << err.what() << std::endl;
        return 1;
    }

    //parse_cmdline(argc, argv, cfg);

    shimmer::infrastructure infra;

    int err = shimmer::load(cfg.database, infra);
    if (err != SHIMMER_SUCCESS) {
        std::cerr << "Problem detected while loading DB" << std::endl;
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
    auto sol_full  = ts1.solution_full();
    std::cout << sol_full << std::endl;

    std::cout << sol_full.rows() << " " << sol_full.cols() << std::endl;

    auto Pbegin = 0;
    auto Plen = num_stations(infra);
    shimmer::save_pressures(cfg.database, infra,
        sol_full.block(0, Pbegin, sol_full.rows(), Plen)
    );

    auto Lbegin = num_stations(infra);
    auto Llen = num_pipes(infra);
    shimmer::save_flowrates(cfg.database, infra,
        sol_full.block(0, Lbegin, sol_full.rows(), Llen)
    );

    return 0;
}