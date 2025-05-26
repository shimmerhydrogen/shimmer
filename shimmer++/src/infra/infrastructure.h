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

#pragma once

#include "sqlite/sqlite.hpp"
#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "solver/conservation_matrices.h"
#include "solver/fluid_solver.h"
#include "solver/time_solver.h"
#include "errors.h"

namespace shimmer {

struct infrastructure {
    infrastructure_graph                    graph;

    optvector<int>                          s_u2i;
    std::vector<int>                        s_i2u;
    optvector<vertex_descriptor>            s_u2vd;
    std::vector<vertex_descriptor>          s_i2vd;

    std::vector<setting_outlet>             settings_outlet;
    std::vector<setting_entry_p_reg>        settings_entry_p_reg;
    std::vector<setting_entry_l_reg>        settings_entry_l_reg;
    std::vector<setting_exit_l_reg>         settings_exit_l_reg;
    std::vector<setting_pipe>               settings_pipe;
    std::vector<setting_compr_stat>         settings_compr_stat;
    std::vector<setting_red_stat>           settings_red_stat;
    std::vector<gas_mass_fractions>         mass_fractions;

    std::vector<station_initial_condition>  sics;
    std::vector<pipe_initial_condition>     pics;
};

int load(const std::string db_filename, infrastructure& infra);

int num_stations(const infrastructure&);
int num_pipes(const infrastructure&);
variable initial_guess(const infrastructure&);

int save_pressures(const std::string&, const infrastructure&, const matrix_t&);
int save_flowrates(const std::string&, const infrastructure&, const matrix_t&);

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

int launch_solver(const config&);

} //namespace shimmer