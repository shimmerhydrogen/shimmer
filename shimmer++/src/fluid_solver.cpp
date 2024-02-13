#include <Eigen/SparseLU>
#include "../src/fluid_solver.h"
#include <iomanip>
namespace shimmer{



double linearized_fluid_solver::temperature(){return Tm_;}
vector_t linearized_fluid_solver::pressure_nodes(){return press_;}
vector_t linearized_fluid_solver::pressure_pipes(){return press_pipes_;}
matrix_t linearized_fluid_solver::x_nodes(){return x_nodes_;}
matrix_t linearized_fluid_solver::x_pipes(){return x_pipes_;}
const incidence& linearized_fluid_solver::get_incidence(){return inc_;}

linearized_fluid_solver::linearized_fluid_solver(
                        double tolerance, 
                        double dt,
                        double Tm,
                        const incidence & inc,
                        const infrastructure_graph& graph):
                        tolerance_(tolerance), dt_(dt), Tm_(Tm),
                        inc_(inc), graph_(graph)  
{
    MAX_ITERS_ = 500;

    num_pipes_ = num_edges(graph_);
    num_nodes_ = num_vertices(graph_);

    /// Underrelaxation coefficients to help convergence
    a_G_ = 0.0;
    a_p_ = 0.0;

    x_nodes_ = build_x_nodes(graph_);
    x_pipes_ = inc_.matrix_in().transpose() * x_nodes_;

    press_.resize(num_nodes_);
    flux_.resize(num_pipes_);
    press_pipes_.resize(num_pipes_);
}



pair_trip_vec_t
linearized_fluid_solver::continuity(
            const vector_t& pressure,
            const vector_t& pressure_old,
            const vector_t& c2)
{
    size_t num_nodes_ = num_vertices(graph_); 
    size_t num_pipes_ = num_edges(graph_);

    vector_t phi_vec = phi_vector(dt_, c2, graph_);
    auto t_sPHI = build_triplets( phi_vec);
    auto t_sA   = build_triplets( inc_.matrix(), 0, num_nodes_ );

    std::vector<triplet_t> triplets =  t_sPHI; 
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());

    vector_t rhs = phi_vec.array() * pressure_old.array();

    return std::make_pair(triplets, rhs);
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

    vector_t rf = resistance_friction(Tm_, c2, flux, graph_);
    vector_t ri = resistance_inertia(dt_, pipes_pressure, inc_, graph_);
    vector_t r = (-rf-ri).array() / ADP_p.array();
    auto t_sR   = build_triplets( r, num_nodes_, num_nodes_);
    auto t_sADP = build_triplets( sADP,  num_nodes_, 0);

    //std::vector<triplet_t> triplets =  t_sAPA;
    std::vector<triplet_t> triplets =  t_sADP;
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());

    vector_t rhs = ((-0.5) * rf.array() * flux.array()
                      - ri.array() * flux_old.array())/ADP_p.array();


    return std::make_pair(triplets, rhs);
}


pair_trip_vec_t
linearized_fluid_solver::boundary(const vector_t& area_pipes,
         const vector_t& p_in,
         const vector_t& flux,
         const vector_t& flux_ext,
         const vector_t& inlet_nodes,
         equation_of_state *eos)
{
    auto rho  = eos->density_correction();

    /// vel [m/s] velocity of the gas within pipes.
    vector_t vel = flux.cwiseQuotient(area_pipes.cwiseProduct(rho));

    sparse_matrix_t sId (num_nodes_, num_nodes_);
    sId.setIdentity();
    auto triplets = build_triplets(sId, 0, num_nodes_ + num_pipes_);

    vector_t rhs = flux_ext;

    for(size_t i = 0; i < inlet_nodes.size(); i++)
    {
        size_t idx = inlet_nodes(i);

        if (vel(idx) > 0.0)
        {
            // Change equation to impose pressure instead of flux
            triplets.push_back(triplet_t(num_pipes_ + num_nodes_ + idx, 0, 1.));
            sId.coeffRef(idx, idx) = 0.0;
            rhs(idx) = p_in(i);
        }
        else
            rhs(idx) = 0.0;

    }

    auto t_sId_bcnd = build_triplets(sId, num_nodes_+num_pipes_, num_nodes_+num_pipes_);
    triplets.insert(triplets.begin(), t_sId_bcnd.begin(), t_sId_bcnd.end());

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


bool linearized_fluid_solver::convergence(const vector_t& diff, 
                const vector_t& sol)
{
    /*  Warning: Marco computes diff = LHS * sol - rhs, with rhs evaluated
        using the current flux (G_k). Here, the rhs is evaluated in the 
        previous step G_{k-1}. 
        res = norm(ADP*p_k -  (Rf * abs(G_k)*G_k) + Ri. * (G_k - G_n)));
            = norm(ADP*p_k -(2 Rf * abs(G_k)+ Ri)*G_k + (Rf *abs(G_k_1)G_k_1 + Ri * G_n))); 
    */
    // diff = lhs sol - rhs

    auto norm_mass = diff.head(num_nodes_).norm();
    auto norm_mom  = diff.segment(num_nodes_, num_pipes_).norm();
    auto residual  = std::max(norm_mass, norm_mom);

    vector_t press_next = a_p_*press_+(1.0-a_p_)*sol.head(num_nodes_);
    vector_t flux_next  = a_G_*flux_ +(1.0-a_G_)*sol.segment(num_nodes_, num_pipes_);

    press_ = press_next;
    flux_  = flux_next;

    std::cout << "sol: ("<< norm_mass<< " , " << norm_mom << ")" <<  std::endl ;
    for (int k = 0; k < sol.size(); ++k)
        std::cout << std::setprecision(16) << sol(k)<< std::endl ;
    std::cout <<  std::endl ;


    if(residual < tolerance_)
        return true; 
    return false;
}



void
linearized_fluid_solver::run(const vector_t& area_pipes,
                            const vector_t& inlet_nodes,
                            const vector_t& p_in,
                            const vector_t& flux_ext,
                            variable& var_time,
                            equation_of_state *eos)
{
    // Initialization of varibales with solution in time n; 
    press_ = var_time.pressure;
    flux_  = var_time.flux;

    std::cout<< "Initializing ..."<< std::endl;

    // Here it takes the x_nodes_/x_pipes_ if needed. For example, for gerg, here it
    // takes these members from lfs, to compute the gerg_params.  
    // Estos se hacen adentro de eos
    //    const gerg_params& gerg_nodes,
    //    const gerg_params& gerg_pipes,


    eos->initialization(this);

    for(size_t iter = 0; iter <= MAX_ITERS_; iter++)
    {   
        std::cout<< "Solver at iteration k ..."<< iter << std::endl;

        press_pipes_ = average(press_, inc_);

        auto [c2_nodes, c2_pipes] = eos->speed_of_sound(this); 

        auto mass = continuity(press_, var_time.pressure, c2_nodes);
        auto mom  = momentum(press_, press_pipes_, flux_, var_time.flux, c2_pipes);       
        auto bcnd = boundary(area_pipes, p_in, flux_,flux_ext, inlet_nodes, eos);
        auto [LHS, rhs]= assemble(mass, mom, bcnd);

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

        vector_t diff = LHS * sol - rhs; 
        if (convergence(diff, sol))
        {
            var_time.pressure = press_;
            var_time.flux = flux_;
            var_time.L_rate = sol.tail(num_nodes_);
            return; 
        }
        
    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;
}





} //end namespace shimmer

