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

#include <Eigen/SparseLU>
#include <iomanip>
#include <fstream>
#include "solver/fluid_solver.h"

namespace shimmer
{

linearized_fluid_solver::linearized_fluid_solver(
                        size_t at_step,
                        const bool& unsteady,
                        double tol,
                        double dt,
                        double Tm,
                        const vector_t& mu,
                        const incidence & inc,
                        const infrastructure_graph& g): at_step_(at_step),
                        is_unsteady_(unsteady),
                        tolerance_(tol), dt_(dt), Tm_(Tm), mu_(mu),
                        inc_(inc), graph_(g)
{
    MAX_ITERS_ = 500;

    num_pipes_ = num_edges(graph_);
    num_nodes_ = num_vertices(graph_);

    /// Underrelaxation coefficients to help convergence
    a_G_ = 0.0;
    a_p_ = 0.0;

    x_nodes_ = build_x_nodes(graph_);
    x_pipes_ = inc_.matrix_in().transpose() * x_nodes_;
}



pair_trip_vec_t
linearized_fluid_solver::continuity(
            const vector_t& pressure_old,
            const vector_t& c2)
{
    size_t num_nodes_ = num_vertices(graph_);
    size_t num_pipes_ = num_edges(graph_);

    vector_t phi_vec = is_unsteady_? phi_vector(dt_, c2, graph_)
                                   : vector_t::Zero(num_nodes_);
    auto t_sPHI = build_triplets( phi_vec);
    auto t_sA   = build_triplets( inc_.matrix(), 0, num_nodes_ );

    std::vector<triplet_t> triplets =  t_sPHI;
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());

    vector_t rhs = phi_vec.array() * pressure_old.array();

    return std::make_pair(triplets, rhs);
}


void
linearized_fluid_solver::impose_edge_station_model(
                        const vector_t& c2_pipes,
                        const vector_t& nodes_pressure,
                        const vector_t& flux,                         
                        sparse_matrix_t& sADP,
                        vector_t& r_scale,
                        vector_t& rhs_mom)

{
    size_t offset = num_nodes_;

    int idx = 0;
    auto edge_range = boost::edges(graph_);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    for (auto itor = begin; itor != end; itor++, idx++)
    {
        auto pipe = graph_[*itor];
        if (pipe.type == pipe_type::PIPE)
            continue;

        const auto& st = pipe.pipe_station;

        auto pipe_num = pipe.branch_num;

        auto s = boost::source(*itor, graph_);
        auto t = boost::target(*itor, graph_);

        auto source_num = graph_[s].i_snum;
        auto target_num = graph_[t].i_snum;

        st->fill_model(st->mode_,
                                      pipe_num,
                                      source_num,
                                      target_num,
                                      var_,
                                      c2_pipes);
    
        // Check if values respect limits, otherwise they are modified
        st->control_hard();

        size_t row = pipe_num + offset;

        //triplets_mom.push_back(triplet_t(row, source_num, st.model_c1()));
        //triplets_mom.push_back(triplet_t(row, target_num, st.model_c2()));

        //WARNING: THIS MUST BE DONE IN ANOTHER WAY!!! TEMPORALLY MODIFYING THE SPARSE MATRIX
        // See Issue #25
        sADP.coeffRef(pipe_num, source_num) = st->model_c1();
        sADP.coeffRef(pipe_num, target_num) = st->model_c2();

        r_scale(pipe_num) = st->model_c3();
        rhs_mom(pipe_num) = st->model_rhs();
    }

    return;
}


pair_trip_vec_t
linearized_fluid_solver::momentum(
        const vector_t& nodes_pressure,
        const vector_t& pipes_pressure,
        const vector_t& flux,
        const vector_t& flux_old,
        const vector_t& c2_nodes,
        const vector_t& c2_pipes)
    {
    size_t num_nodes_ = num_vertices(graph_);
    size_t num_pipes_ = num_edges(graph_);

    sparse_matrix_t sADP = adp_matrix(c2_pipes, graph_, inc_);
    vector_t ADP_p = sADP.cwiseAbs() * nodes_pressure;

    vector_t rf = resistance_friction(Tm_, mu_, c2_pipes, flux, graph_);
    vector_t r = -rf;
    vector_t rhs = -0.5 * rf.cwiseProduct(flux);

    if (is_unsteady_)
    {
        vector_t ri = resistance_inertia(dt_, pipes_pressure, inc_, graph_);
        r -= ri;
        rhs -= ri.cwiseProduct(flux_old);
    }

    vector_t r_scale = r.cwiseQuotient(ADP_p);
    vector_t rhs_scale = rhs.cwiseQuotient(ADP_p);

    impose_edge_station_model(c2_pipes, nodes_pressure, flux,
                              sADP,r_scale, rhs_scale);

    auto t_sR   = build_triplets( r_scale , num_nodes_, num_nodes_);
    auto t_sADP = build_triplets( sADP,  num_nodes_, 0);

    std::vector<triplet_t> triplets =  t_sADP;
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());

    return std::make_pair(triplets, rhs_scale);
}


pair_trip_vec_t
linearized_fluid_solver::boundary(const vector_t& area_pipes,
         const vector_t& flux,
         equation_of_state *eos)
{

    size_t offset =  num_nodes_ + num_pipes_;

    // This is here due to the first code in matlab. But could be not necessary.
    auto rho  = eos->density();
    /// vel [m/s] velocity of the gas within pipes.
    vector_t vel = flux.cwiseQuotient(area_pipes.cwiseProduct(rho));


    sparse_matrix_t sId (num_nodes_, num_nodes_);
    sId.setIdentity();
    auto triplets = build_triplets(sId, 0, offset);

    vector_t rhs(num_nodes_);

    int idx = 0;
    auto v_range = boost::vertices(graph_);
    for(auto itor = v_range.first; itor != v_range.second; itor++, idx++)
    {
        if (not graph_[*itor].node_station) {
            std::string errstr = "Graph node " + std::to_string(*itor) +
                " points to an invalid station"; 
            throw std::invalid_argument(errstr);
        }
        auto bnd =  graph_[*itor].node_station->boundary();

        switch(bnd.type())
        {
            case constraint_type::L_EQUAL:
                triplets.push_back(triplet_t(idx + offset, idx + offset, 1.));
                break;
            case constraint_type::P_EQUAL:
                triplets.push_back(triplet_t(idx + offset, idx, 1.));
                break;
            default:
                throw std::invalid_argument("Boundary type not valid.");
        }
        rhs(idx) = bnd.value(at_step_);
    }

    return  std::make_pair(triplets, rhs);
}


/*
    num_nodes        | PHI    A   Ic| | p |
    num_pipes        | APA   -R   0 | | G |
    num_nodes_bnd    |  Ib    0   Ib| | Gb|
*/


std::pair<sparse_matrix_t, vector_t>
linearized_fluid_solver::assemble(
        const pair_trip_vec_t& lhs_rhs_mass,
        const pair_trip_vec_t& lhs_rhs_mom,
        const pair_trip_vec_t& lhs_rhs_bnd)
{
    std::vector<triplet_t> triplets = lhs_rhs_mass.first;
    auto lhs_mom_begin = lhs_rhs_mom.first.cbegin();
    auto lhs_mom_end   = lhs_rhs_mom.first.cend();
    auto lhs_bnd_begin = lhs_rhs_bnd.first.cbegin();
    auto lhs_bnd_end   = lhs_rhs_bnd.first.cend();

    triplets.insert(triplets.begin(),lhs_mom_begin, lhs_mom_end);
    triplets.insert(triplets.begin(),lhs_bnd_begin, lhs_bnd_end);

    size_t system_size = num_nodes_ + num_pipes_ + num_nodes_;

    sparse_matrix_t LHS(system_size, system_size);
    LHS.setFromTriplets(triplets.begin(), triplets.end());

    vector_t rhs = vector_t::Zero(system_size);
    rhs.head(num_nodes_) =  lhs_rhs_mass.second;
    rhs.segment(num_nodes_, num_pipes_) = lhs_rhs_mom.second;
    rhs.tail(num_nodes_) =  lhs_rhs_bnd.second;
    return std::make_pair(LHS, rhs);
}


bool linearized_fluid_solver::convergence(
                const vector_t& sol)
{
    vector_t diff = sol - var_.make_vector();
    auto norm_mass = diff.head(num_nodes_).norm()/(var_.pressure).norm();
    auto norm_mom  = diff.segment(num_nodes_, num_pipes_).norm()/(var_.flux).norm();
    auto residual  = std::max(norm_mass, norm_mom);

    //vector_t press_next = a_p_*var_.pressure+(1.0-a_p_)*sol.head(num_nodes_);
    //vector_t flux_next  = a_G_*var_.flux +(1.0-a_G_)*sol.segment(num_nodes_, num_pipes_);

    var_.pressure = sol.head(num_nodes_);//press_next;
    var_.flux  = sol.segment(num_nodes_, num_pipes_);//flux_next;
    var_.L_rate = sol.tail(num_nodes_);

    if(residual < tolerance_)
        return true;
    return false;
}


bool
linearized_fluid_solver::run(const vector_t& area_pipes,
                            const variable& var_guess,
                            const variable& var_time,
                            equation_of_state *eos,                        
                            size_t at_iteration)
{
    press_pipes_.resize(num_pipes_);

    // Initialization of variables with solution in time n;
    var_.pressure = var_guess.pressure;
    var_.flux  = var_guess.flux;
    var_.L_rate = vector_t::Zero(num_nodes_);//flux_ext;
    std::cout<< "Initializing ..."<< std::endl;

    eos->initialization(this);

    for(size_t iter = 0; iter <= MAX_ITERS_; iter++)
    {
        std::cout << "---------------------------------------------------"<< std::endl;
        std::cout<< "Fluid solver at iteration k ..............."<< iter << std::endl;
        std::cout << "---------------------------------------------------"<< std::endl;


        press_pipes_ = average(var_.pressure, inc_);
        auto [c2_nodes, c2_pipes] = eos->speed_of_sound(this);

        auto mass = continuity(var_time.pressure, c2_nodes);
        auto mom  = momentum(var_.pressure, press_pipes_, var_.flux, var_time.flux, c2_nodes, c2_pipes);
        auto bcnd = boundary(area_pipes, var_.flux, eos);
        auto [LHS, rhs] = assemble(mass, mom, bcnd);


        Eigen::SparseLU<sparse_matrix_t> solver;
        solver.compute(LHS);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error factorizing LHS" <<std::endl;
#if 0
            size_t count = 0;
            for (int k = 0; k < LHS.outerSize(); ++k)
            {
                for (itor_t it(LHS,k); it; ++it, count++)
                {
                    std::cout << std::setprecision(16) << "" << it.row()
                                << " , " << it.col() << " , " << it.value()
                                << " ; " << std::endl ;
                }
            }
#endif
            exit(1);
        }

        vector_t sol = solver.solve(rhs);


        if(solver.info() != Eigen::Success) {
            std::cout << "Error solving system" <<std::endl;
            exit(1);
        }

        std::string str_iter =  std::to_string(at_iteration); 
        std::string str_step =  std::to_string(at_step_); 
        
        if (convergence(sol))
        {
            c2_nodes_ = c2_nodes;
            c2_pipes_ = c2_pipes;

            return true;
        }

    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;
    return false;
}


bool
linearized_fluid_solver::check_hard_constraints(size_t step)
{
    bool pass_all = true;
    int i = 0;
    auto v_range = boost::vertices(graph_);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
    {
        bool pass = graph_[*itor].node_station->check_hard(var_.pressure[i], var_.L_rate[i], step);
        std::cout<< " * Hard ("<< i << ") : "<< pass << std::endl;
        pass_all = pass_all && pass;
    }

    return pass_all;
}


bool
linearized_fluid_solver::check_hard_controls(size_t step)
{
    bool pass_all = true;

    auto edge_range = boost::edges(graph_);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    for (auto itor = begin; itor != end; itor++)
    {
        const auto& pipe = graph_[*itor];

        if (pipe.type == pipe_type::PIPE)
            continue;

        auto& st = pipe.pipe_station;

        if (!st->is_on()) continue;

        int pipe_num = pipe.branch_num;

        auto s = boost::source(*itor, graph_);
        auto t = boost::target(*itor, graph_);

        int source_num = graph_[s].i_snum;
        int target_num = graph_[t].i_snum;

        // =================================================================
        // fill all controls for verification:
        for (edge_station::control::mode& m : st->controls_on)
            st->fill_model(m, pipe_num, source_num, target_num, var_, c2_pipes_);
        

        // =================================================================

        /* Check other controls mode:

        Once the modes are filled with the current values, the constraints
        are checked.
        If one mode constraint is not verified then run simulation again using
        the new control mode
        */
        size_t idx = 0;
        bool pass;

        for (const auto& m : st->controls_on)
        {
            pass = m.check_hard();
            if (!pass)
            {
                st->change_mode_on(idx);
                break;
            }
            idx++;
        }
        // =================================================================
        std::cout << " * Hard (" << pipe_num << ") : " << pass << std::endl;

        pass_all = pass_all && pass;
    }

    return pass_all;
}


void
linearized_fluid_solver::check_soft_constraints(size_t step)
{
    int i = 0;
    auto v_range = boost::vertices(graph_);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        graph_[*itor].node_station->check_soft(var_.pressure[i], var_.L_rate[i], step);
}

/*
void
linearized_fluid_solver::check_soft_controls(size_t step)
{
    auto edge_range = boost::edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    for (auto itor = begin; itor != end; itor++)
    {
        auto pipe = g[*itor];

        if (pipe.type == edge_type::pipe) continue;

        graph_[*itor].node_station->check_soft(var_.pressure[i], var_.L_rate[i], step);
    }
}
*/

bool
linearized_fluid_solver::check_constraints(size_t step)
{
    check_soft_constraints(step);

    return check_hard_constraints(step);
}


bool
linearized_fluid_solver::check_controls(size_t step)
{
    if (check_hard_controls(step))
    {
        //check_soft_controls(step);
        return true;

    }

    return false;
}


} //end namespace shimmer
