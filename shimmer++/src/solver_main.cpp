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
#include <getopt.h>

#include "sol/sol.hpp"


#include "infra/infrastructure.h"
#include "errors.h"
#include "infra/launch_solver.h"


static void
register_usertypes(sol::state& lua)
{
    using namespace shimmer;
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
    config_ut["refine"] = &config::refine;
    config_ut["dx"] = &config::dx;
    config_ut["quality_tracking"] = &config::quality_tracking;
    config_ut["qt_steady"] = &config::qt_steady;
}

int main(int argc, char **argv)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);


    int c;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"init-db",     required_argument, 0,  'i'},
            {0,             0,                 0,   0 }
        };

        c = getopt_long(argc, argv, "i:",
                        long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            //case 0:
            case 'i':
                std::cout << "Initializing database file '";
                std::cout << optarg << "'" << std::endl;
                return shimmer::initdb(optarg);
                break;

            default:
                std::cout << "Unknown option specified" << std::endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;

    if (argc < 1) {
        std::cerr << "Please specify configuration file name" << std::endl;
        return 1;
    }

    sol::state lua;
    lua.open_libraries(
        sol::lib::base, sol::lib::math, sol::lib::io,
        sol::lib::table, sol::lib::string
        );

    register_usertypes(lua);    

    shimmer::config cfg;
    cfg.refine = false;
    cfg.dx = 100e3;
    cfg.quality_tracking = false;
    cfg.qt_steady = false;
    lua["config"] = &cfg;

    auto sresult = lua.safe_script_file(argv[0], sol::script_pass_on_error);
    if (!sresult.valid()) {
        sol::error err = sresult;
        std::cerr << "Solver configuration: " << err.what() << std::endl;
        return 1;
    }

    if (cfg.qt_steady)
        cfg.steps = 1; 

    if (cfg.quality_tracking)
        return launch_solver_qt(cfg);
    else
        return launch_solver(cfg);
}