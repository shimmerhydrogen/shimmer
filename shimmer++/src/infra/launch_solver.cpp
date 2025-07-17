
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
    }

    assert( num_stations(infra) != 0 );
    assert( num_pipes(infra) != 0 );

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

    // Solver
    using time_solver_t = shimmer::time_solver<shimmer::papay,
        shimmer::viscosity_type::Constant>;

    time_solver_t ts1(infra.graph, cfg.temperature);
    ts1.initialization(guess, cfg.dt_std, cfg.tol_std, y_nodes, y_pipes);  
    ts1.advance(cfg.dt, cfg.steps, cfg.tol, y_nodes, y_pipes);
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
    }

    assert( num_stations(infra) != 0 );
    assert( num_pipes(infra) != 0 );

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

    // Solver
    using solver_t = shimmer::qt_solver<shimmer::papay,
        shimmer::viscosity_type::Constant>;

    solver_t qt(infra, cfg.temperature);
    qt.initialization(guess, cfg.dt_std, cfg.tol_std, y_nodes, y_pipes);  
    //qt.advance(cfg.dt, cfg.steps, cfg.tol, y_nodes, y_pipes);
    auto sol_full  = qt.solution_full();
    auto vel_full  = qt.velocity_full();

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

} //end namespace shimmer