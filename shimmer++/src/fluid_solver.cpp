#include <Eigen/SparseLU>
#include "../src/fluid_solver.h"
#include "../src/nonpipe_pipe.h"
#include <iomanip>
#include <fstream>

namespace shimmer
{



double linearized_fluid_solver::temperature(){return Tm_;}
vector_t linearized_fluid_solver::pressure_nodes(){return var_.pressure;}
vector_t linearized_fluid_solver::pressure_pipes(){return press_pipes_;}
matrix_t linearized_fluid_solver::x_nodes(){return x_nodes_;}
matrix_t linearized_fluid_solver::x_pipes(){return x_pipes_;}
const incidence& linearized_fluid_solver::get_incidence(){return inc_;}

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
    /*
    std::cout << "* PHI: " <<  std::endl;
    for (int k = 0; k < phi_vec.size(); ++k)
        std::cout << std::setprecision(16) <<phi_vec(k)<< std::endl;
    std::cout <<  std::endl ;


    std::ofstream ofs;
    ofs.open("p_n_karol.dat", std::ios_base::app);

    for(size_t j = 0; j < pressure_old.size(); j++)
        ofs << std::setprecision(16) << pressure_old(j) << " ";
    ofs << std::endl;
    ofs.close();

    std::cout << "* p_n: " <<  std::endl;
    for (int k = 0; k < pressure_old.size(); ++k)
        std::cout << std::setprecision(16) << pressure_old(k)<< std::endl;
    std::cout <<  std::endl;
    */

    return std::make_pair(triplets, rhs);
}


void
linearized_fluid_solver::impose_edge_station_model(std::vector<triplet_t>& triplets_mom,
    vector_t& rhs_mom)
{
    size_t offset = num_nodes_;

    int idx = 0;
    auto edge_range = boost::edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    for (auto itor = begin; itor != end; itor++, idx++)
    {
        auto pipe = g[*itor];
        if (pipe.type == edge_type::pipe) continue;

        auto& st = pipe.pipe_station;

        //if (!st->active(at_step_))  continue;

        auto pipe_idx = pipe.branch_num;

        size_t row = pipe_idx + offset;
        auto s = source(*itor, g);
        auto t = target(*itor, g);

        auto source_node = g[s].node_num;
        auto target_node = g[t].node_num;

        triplets_mom.push_back(triplet_t(row, source_node, st.model_c1()));
        triplets_mom.push_back(triplet_t(row, target_node, st.model_c2()));
        triplets_mom.push_back(triplet_t(row, row, st.model_c3()));
        rhs_mom(row) = st.model_rhs();
    }

    return;
}


void
linearized_fluid_solver::control_stations(
            const vector_t& c2_nodes,
            const vector_t& nodes_pressure,
            const vector_t& pipes_pressure,
            const vector_t& flux)
{
    size_t offset = num_nodes_;

    int idx = 0;
    auto edge_range = boost::edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    for (auto itor = begin; itor != end; itor++, idx++)
        {
        auto pipe = g[*itor];
        if (pipe.type == edge_type::pipe) continue;

        auto& st = pipe.pipe_station;

        auto pipe_idx = pipe.branch_num;

        size_t row = pipe_idx + offset;
        auto s = source(*itor, g);
        auto t = target(*itor, g);

        auto source_node = g[s].node_num;
        auto target_node = g[t].node_num;

        //for (auto c = 0; c < st.num_controls();c++)
        //{
        auto p_in = nodes_pressure[source_node];
        auto p_out = nodes_pressure[target_node];

        // Set variable
        switch (st->which_control_type(at_step_))
        {
            case control::type::SHUT_OFF:
            case control::type::BY_PASS:
                set.set_c3(p_in < p_out);
                break;
            case control::type::BETA:
                auto beta = p_out /p_in;
                st.set_c1(beta);
                break;
            case control::type::POWER_DRIVER:
            {
                auto gamma = 1.4; // Or read from GERG
                auto ck = gamma - 1.0 / gamma;
                auto beta = p_out /p_in;
                //beta = edge_stations::compressor_beta(p_in, p_out, st.internals());
                auto ZTR = c2_nodes[source_node];
                auto K = ZTR / st.efficiency();
                auto G = flux[pipe_idx];
                auto KGB = K * G * beta;

                auto c3 = (K / ck) * (std::pow(beta, ck) - 1.0);
                auto pwd = c3 * G;

                st.set_c1(-KGB / p_in);
                st.set_c2( KGB / p_out);
                st.set_c3(c3);
                st.set_rhs(pwd);
                break;
            }
            case control::type::PRESSURE_IN:
                st.set_rhs(p_in);
                break;
            case control::type::PRESSURE_OUT:
                st.set_rhs(p_out);
                break;
            case control::type::FLUX:
                st.set_rhs(flux[pipe_idx]);
                break;
            default:
                std::cout << "ERROR: Fluid solver does not know this control type.\n";
                throw std::exception;
        }

        // Control variable
        st.control_hard( at_step_* dt_);

//      }
    }

    return;
}

pair_trip_vec_t
linearized_fluid_solver::momentum(
         const vector_t& nodes_pressure,
         const vector_t& pipes_pressure,
         const vector_t& flux,
         const vector_t& flux_old,
         const vector_t& c2)
{
    size_t num_nodes_ = num_vertices(graph_);
    size_t num_pipes_ = num_edges(graph_);
    //size_t num_pipes_ext = num_pipes_;

    sparse_matrix_t sADP = adp_matrix(c2, graph_, inc_);
    vector_t ADP_p = sADP.cwiseAbs() * nodes_pressure;

    vector_t rf = resistance_friction(Tm_, mu_, c2, flux, graph_);
    vector_t r = -rf;
    vector_t rhs = -0.5 * rf.cwiseProduct(flux);

    if (is_unsteady_)
    {
        vector_t ri = resistance_inertia(dt_, pipes_pressure, inc_, graph_);
        r -= ri;
        rhs -= ri.cwiseProduct(flux_old);

        /*
        std::cout << "* ri: " <<  std::endl;
        for (int k = 0; k < ri.size(); ++k)
            std::cout << std::setprecision(16) <<ri(k)<< std::endl;
        std::cout <<  std::endl ;
        */
    }

    vector_t r_scale = r.cwiseQuotient(ADP_p);
    auto t_sR   = build_triplets( r_scale , num_nodes_, num_nodes_);
    auto t_sADP = build_triplets( sADP,  num_nodes_, 0);

    std::vector<triplet_t> triplets =  t_sADP;
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());

    vector_t rhs_scale = rhs.cwiseQuotient(ADP_p);

    impose_edge_station_model(triplets, rhs_scale);
    /*
    std::cout << "* nodes_pressure: " <<  std::endl;
    for (int k = 0; k < nodes_pressure.size(); ++k)
        std::cout << std::setprecision(16) <<nodes_pressure(k)<< std::endl;
    std::cout <<  std::endl ;

    std::cout << "* rf: " <<  std::endl;
    for (int k = 0; k < rf.size(); ++k)
        std::cout << std::setprecision(16) <<rf(k)<< std::endl;
    std::cout <<  std::endl ;


    std::cout << "* c2_pipes: " <<  std::endl;
    for (int k = 0; k < c2.size(); ++k)
        std::cout << std::setprecision(16) <<c2(k)<< std::endl;
    std::cout <<  std::endl ;


    std::cout << "* ADP_p: " <<  std::endl;
    for (int k = 0; k < ADP_p.size(); ++k)
        std::cout << std::setprecision(16) <<ADP_p(k)<< std::endl;
    std::cout <<  std::endl ;

    std::cout << "* r_scale: " <<  std::endl;
    for (int k = 0; k < r_scale.size(); ++k)
        std::cout << std::setprecision(16) <<r_scale(k)<< std::endl;
    std::cout <<  std::endl ;

    std::cout << "* rhs_scale: " <<  std::endl;
    for (int k = 0; k < rhs_scale.size(); ++k)
        std::cout << std::setprecision(16) <<rhs_scale(k)<< std::endl;
    std::cout <<  std::endl ;
    */

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

/*
bool linearized_fluid_solver::convergence_denerg(
                const vector_t& sol)
{

}
*/

bool linearized_fluid_solver::convergence(
                const vector_t& sol)
{
    /*  Warning: Marco computes diff = LHS^(k+1) * sol^(k+1) - rhs^(k+1)
        While here, sol^{k+1} - sol^{k}
    */

    vector_t diff = sol - var_.make_vector();
    auto norm_mass = diff.head(num_nodes_).norm()/(var_.pressure).norm();
    auto norm_mom  = diff.segment(num_nodes_, num_pipes_).norm()/(var_.flux).norm();
    auto residual  = std::max(norm_mass, norm_mom);

    vector_t press_next = a_p_*var_.pressure+(1.0-a_p_)*sol.head(num_nodes_);
    vector_t flux_next  = a_G_*var_.flux +(1.0-a_G_)*sol.segment(num_nodes_, num_pipes_);

    var_.pressure = press_next;
    var_.flux  = flux_next;
    var_.L_rate = sol.tail(num_nodes_);
/*
    std::cout << "sol: ("<< norm_mass<< " , " << norm_mom << ")" <<  std::endl ;
    for (int k = 0; k < sol.size(); ++k)
        std::cout << std::setprecision(16) << sol(k)<< std::endl ;
    std::cout <<  std::endl ;
*/

    if(residual < tolerance_)
        return true;
    return false;
}



bool
linearized_fluid_solver::run(const vector_t& area_pipes,
                            //const vector_t& flux_ext,
                            const variable& var_guess,
                            const variable& var_time,
                            equation_of_state *eos)
{
    press_pipes_.resize(num_pipes_);

    // Initialization of varibales with solution in time n;
    var_.pressure = var_guess.pressure;
    var_.flux  = var_guess.flux;
    var_.L_rate = vector_t::Zero(num_nodes_);//flux_ext;
    std::cout<< "Initializing ..."<< std::endl;

    // Here it takes the x_nodes_/x_pipes_ if needed. For example, for gerg, here it
    // takes these members from lfs, to compute the gerg_params.
    // Estos se hacen adentro de eos
    //    const gerg_params& gerg_nodes,
    //    const gerg_params& gerg_pipes,


    eos->initialization(this);

    for(size_t iter = 0; iter <= MAX_ITERS_; iter++)
    {
        std::cout<< "---------------------------------" << std::endl;
        std::cout<< "Solver at iteration k ..."<< iter << std::endl;

        press_pipes_ = average(var_.pressure, inc_);
        auto [c2_nodes, c2_pipes] = eos->speed_of_sound(this);

        control_edge_stations(c2_nodes);

        auto mass = continuity(var_time.pressure, c2_nodes);
        auto mom  = momentum(var_.pressure, press_pipes_, var_.flux, var_time.flux, c2_pipes);
        auto bcnd = boundary(area_pipes, var_.flux, eos);
        auto [LHS, rhs] = assemble(mass, mom, bcnd);

        /*
        //std::vector<triplet_t> LHS_MOM  = mom.second;
        std::cout << "Momemtum " << std::endl;
        for (const triplet_t & trip : mom.first)
        {
            std::cout << std::setprecision(16) << "" << trip.row()
                            << " , " << trip.col() << " , " << trip.value()
                            << " ; " << std::endl ;
        }
        */

        Eigen::SparseLU<sparse_matrix_t> solver;
        solver.compute(LHS);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error factorizing LHS" <<std::endl;

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
            exit(1);
        }



        vector_t sol = solver.solve(rhs);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error solving system" <<std::endl;
            exit(1);
        }

        /*
        std::cout << "LHS : " <<std::endl;
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


        std::cout << "rhs = "<< std::endl;
        for (int k = 0; k < rhs.size(); ++k)
            std::cout << "  " << rhs[k]  <<  std::endl;

        std::cout << " * XXX_k at iter ...."<< iter << std::endl;
        for (int k = 0; k < sol.size(); ++k)
            std::cout << "  " << sol[k]  <<  std::endl;
        */
        std::cout<< "Solver at iteration k ..."<< iter << std::endl;

        if (convergence(sol))
        {
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
linearized_fluid_solver::check_edge_hard_constraints(size_t step)
    {}

/*
bool
linearized_fluid_solver::check_edge_control_constraints(size_t step)
{
    bool pass_all = true;

    auto edge_range = boost::edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    for (auto itor = begin; itor != end; itor++)
    {
        auto pipe = g[*itor];

        if (pipe.type == edge_type::pipe) continue;

        auto& st = pipe.pipe_station;

        if (!st->active(step)) continue;

        auto pipe_idx = pipe.branch_num;
        auto s_node_idx = g[source(*itor, g)].node_num;
        auto t_node_idx = g[target(*itor, g)].node_num;

        double p_in  = var_.pressure[s_node_idx];
        double p_out = var_.pressure[t_node_idx];
        double flow  = var_.flux[pipe_idx];
        bool pin_greater_pout = p_in < p_out;

        using cvar_t = control::variable_to_control;

        auto control_vars = std::unordered_map<double>{ {cvar_t::PRESSURE_IN, pin },
                                                        {cvar_t::PRESSURE_OUT, p_out},
                                                        {cvar_t::FLOW, flow},
                                                        {cvar_t::PIN_LOWER_POUT, pin < pout} };

        auto pass = pipe.edge_station->control_hard(control_vars, step * dt_);

        std::cout << " * Hard (" << idx << ") : " << pass << std::endl;

        pass_all = pass_all && pass;
    }

    return pass_all;
}
*/

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
linearized_fluid_solver::check_edge_soft_constraints(size_t step)
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
    if(check_hard_constraints(step))
    {
        check_soft_constraints(step);
        return true;
    }
    return false;
}

/*
bool
linearized_fluid_solver::check_edge_constraints(size_t step)
{

    if(check_edge_hard_constraints(step))
        check_edge_soft_constraints(step);
        return true;
    }
    return false;
}
*/

} //end namespace shimmer
