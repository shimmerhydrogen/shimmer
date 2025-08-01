
#include <filesystem>

#include "infra/launch_solver.h"
#include "errors.h"


namespace fs = std::filesystem;

namespace shimmer {

int launch_solver(const config& cfg)
{
    shimmer::infrastructure infra;
    
    if(cfg.refine) {
        shimmer::infrastructure infrain;
        int err = shimmer::load(cfg.database, infrain);
        if (err != SHIMMER_SUCCESS) {
            std::cerr << "Problem detected while loading DB" << std::endl;
            return 1;
        }
        refine_pipes(infrain, infra, cfg.dx);
    } else {
        int err = shimmer::load(cfg.database, infra);
        if (err != SHIMMER_SUCCESS) {
            std::cerr << "Problem detected while loading DB" << std::endl;
            return 1;
        }
        infra.num_original_stations = num_stations(infra);
        infra.num_original_pipes = num_pipes(infra);

        auto edge_range = boost::edges(infra.graph);
        for (auto itor = edge_range.first; itor != edge_range.second; itor++){
            infra.p_i2ed.push_back(*itor);
        }
    }

    assert( num_stations(infra) != 0 );
    assert( num_pipes(infra) != 0 );

    shimmer::variable guess = initial_guess(infra);

    /* END GAS MASS FRACTIONS */

    // Solver
    using time_solver_t = shimmer::time_solver<shimmer::papay,
        shimmer::viscosity_type::Constant>;

    time_solver_t ts1(infra.graph, cfg.temperature);
    ts1.initialization(guess, cfg.dt_std, cfg.tol_std);  
    ts1.advance(cfg.dt, cfg.steps, cfg.tol);
    auto sol_full  = ts1.solution_full();
    auto vel_full  = ts1.velocity_full();

    // Post-processing
    std::cout << sol_full << std::endl;
    std::cout << sol_full.rows() << " " << sol_full.cols() << std::endl;

    std::string outfile = cfg.database;
    if (cfg.refine) {
        fs::path path(cfg.database);
        std::string dir = path.parent_path();
        std::string file = path.filename();
        outfile =  "refined_" + file;
        if (initdb(outfile) != 0) {
            std::cerr << "Problem creating output db" << std::endl;
            return -1;
        }

        shimmer::store(outfile, infra);
    }

    auto Pbegin = 0;
    auto Plen = num_stations(infra);
    shimmer::save_pressures(outfile, infra,
        sol_full.block(0, Pbegin, sol_full.rows(), Plen)
    );

    auto Gbegin = num_stations(infra);
    auto Glen = num_pipes(infra);
    shimmer::save_flowrates(outfile, infra,
        sol_full.block(0, Gbegin, sol_full.rows(), Glen)
    );

    auto Lbegin = Gbegin + Glen;
    auto Llen = num_stations(infra);
    shimmer::save_flowrates_stations(outfile, infra,
        sol_full.block(0, Lbegin, sol_full.rows(), Llen)
    );

    shimmer::save_velocities(outfile, infra, vel_full);

    return 0;
}

int launch_solver_qt(const config& cfg)
{
    shimmer::infrastructure infra;
    
    if(cfg.refine) {
        shimmer::infrastructure infrain;
        int err = shimmer::load(cfg.database, infrain);
        if (err != SHIMMER_SUCCESS) {
            std::cerr << "Problem detected while loading DB" << std::endl;
            return 1;
        }
        refine_pipes(infrain, infra, cfg.dx);
    } else {
        int err = shimmer::load(cfg.database, infra);
        if (err != SHIMMER_SUCCESS) {
            std::cerr << "Problem detected while loading DB" << std::endl;
            return 1;
        }
        infra.num_original_stations = num_stations(infra);
        infra.num_original_pipes = num_pipes(infra);

        auto edge_range = boost::edges(infra.graph);
        for (auto itor = edge_range.first; itor != edge_range.second; itor++){
            edge_descriptor ed = *itor;
            infra.p_i2ed.push_back(ed);
        }

    }

    assert( num_stations(infra) != 0 );
    assert( num_pipes(infra) != 0 );

    shimmer::variable guess = initial_guess(infra);

    // Solver
    using solver_t = shimmer::qt_solver<shimmer::papay,
        shimmer::viscosity_type::Constant>;

    solver_t qt(infra, cfg.temperature);
    qt.initialization(guess, cfg.dt_std, cfg.tol_std);  
    qt.advance(cfg.dt, cfg.steps, cfg.tol);
    auto sol_full  = qt.solution_full();
    auto vel_full  = qt.velocity_full();
    auto x_full  = qt.molar_fractions_full();

    // Post-processing
    std::cout << sol_full << std::endl;
    std::cout << sol_full.rows() << " " << sol_full.cols() << std::endl;

    std::string outfile = cfg.database;
    if (cfg.refine) {
        fs::path path(cfg.database);
        std::string dir = path.parent_path();
        std::string file = path.filename();
        outfile =  "refined_" + file;
        if (initdb(outfile) != 0) {
            std::cerr << "Problem creating output db" << std::endl;
            return -1;
        }

        shimmer::store(outfile, infra);
    }

    auto Pbegin = 0;
    auto Plen = num_stations(infra);
    shimmer::save_pressures(outfile, infra,
        sol_full.block(0, Pbegin, sol_full.rows(), Plen)
    );

    auto Gbegin = num_stations(infra);
    auto Glen = num_pipes(infra);
    shimmer::save_flowrates(outfile, infra,
        sol_full.block(0, Gbegin, sol_full.rows(), Glen)
    );

    auto Lbegin = Gbegin + Glen;
    auto Llen = num_stations(infra);
    shimmer::save_flowrates_stations(outfile, infra,
        sol_full.block(0, Lbegin, sol_full.rows(), Llen)
    );

    shimmer::save_velocities(outfile, infra, vel_full);

    shimmer::save_molar_fractions(outfile, infra, x_full);

    return 0;
}

} //end namespace shimmer