#include <Eigen/SparseLU>
#include "../src/fluid_solver.h"
#include <iomanip>
namespace shimmer{



double linearized_fluid_solver::temperature(){return Tm_;}
vector_t linearized_fluid_solver::pressure_nodes(){return var_.pressure;}
vector_t linearized_fluid_solver::pressure_pipes(){return press_pipes_;}
matrix_t linearized_fluid_solver::x_nodes(){return x_nodes_;}
matrix_t linearized_fluid_solver::x_pipes(){return x_pipes_;}
const incidence& linearized_fluid_solver::get_incidence(){return inc_;}

linearized_fluid_solver::linearized_fluid_solver(
                        const bool& unsteady,
                        double tol, 
                        double dt,
                        double Tm,
                        const vector_t& mu,
                        const incidence & inc,
                        const infrastructure_graph& g):
                        tolerance_(tol), dt_(dt), Tm_(Tm), mu_(mu),
                        inc_(inc), graph_(g), is_unsteady_(unsteady)  
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

    
    std::cout << "* p_n: " <<  std::endl;
    for (int k = 0; k < pressure_old.size(); ++k)
        std::cout << std::setprecision(16) << pressure_old(k)<< std::endl;
    std::cout <<  std::endl ;
*/
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

    vector_t rf = resistance_friction(Tm_, mu_, c2, flux, graph_);
    vector_t r = -rf; 
    vector_t rhs = -0.5 * rf.cwiseProduct(flux);

    if (is_unsteady_)
    {
        vector_t ri = resistance_inertia(dt_, pipes_pressure, inc_, graph_);
        r -= ri; 
        rhs -= ri.cwiseProduct(flux_old);
    }
    
    vector_t r_scale = r.cwiseQuotient(ADP_p);
    auto t_sR   = build_triplets( r_scale , num_nodes_, num_nodes_);
    auto t_sADP = build_triplets( sADP,  num_nodes_, 0);

    std::vector<triplet_t> triplets =  t_sADP;
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());

    vector_t rhs_scale = rhs.cwiseQuotient(ADP_p);
    return std::make_pair(triplets, rhs_scale);
}


pair_trip_vec_t
linearized_fluid_solver::boundary(const vector_t& area_pipes,
         double p_in,
         const vector_t& flux,
         const vector_t& flux_ext,
         const vector_t& inlet_nodes,
         equation_of_state *eos)
{
    auto rho  = eos->density();

    /// vel [m/s] velocity of the gas within pipes.
    vector_t vel = flux.cwiseQuotient(area_pipes.cwiseProduct(rho));

    sparse_matrix_t sId (num_nodes_, num_nodes_);
    sId.setIdentity();
    auto triplets = build_triplets(sId, 0, num_nodes_ + num_pipes_);

    vector_t rhs = flux_ext;

/*
    std::cout << "flux_ext: " <<  std::endl ;
    for (int k = 0; k < flux_ext.size(); ++k)
        std::cout <<"    " << flux_ext(k)<< std::endl ;
    std::cout <<  std::endl ;
*/

    for(size_t i = 0; i < inlet_nodes.size(); i++)
    {
        size_t idx = inlet_nodes(i);

        /*
        if (vel(idx) > 0.0)
        {
            // Change equation to impose pressure instead of flux

            // check if this is : 
            //                 triplet_t(num_pipes_ + num_nodes_ + idx, idx, 1.))
            triplets.push_back(triplet_t(num_pipes_ + num_nodes_ + idx, 0, 1.));
            sId.coeffRef(idx, idx) = 0.0;
            rhs(idx) = p_in(i);
        }
        else
            rhs(idx) = 0.0;
        */

            // Change equation to impose pressure instead of flux
            triplets.push_back(triplet_t(num_pipes_ + num_nodes_ + idx, idx, 1.));
            sId.coeffRef(idx, idx) = 0.0;
            rhs(idx) = p_in;
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



void
linearized_fluid_solver::run(const vector_t& area_pipes,
                            const vector_t& inlet_nodes,
                            double p_in,
                            const vector_t& flux_ext,
                            const variable& var_guess,
                            variable& var_time,
                            equation_of_state *eos)
{
    press_pipes_.resize(num_pipes_);

    // Initialization of varibales with solution in time n; 
    var_.pressure = var_guess.pressure;
    var_.flux  = var_guess.flux;
    var_.L_rate = flux_ext;
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
        auto mass = continuity(var_time.pressure, c2_nodes);
        auto mom  = momentum(var_.pressure, press_pipes_, var_.flux, var_time.flux, c2_pipes);       
        auto bcnd = boundary(area_pipes, p_in, var_.flux,flux_ext, inlet_nodes, eos);
        auto [LHS, rhs]= assemble(mass, mom, bcnd);

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

/*               
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
*/        

        vector_t sol = solver.solve(rhs);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error solving system" <<std::endl;
            exit(1);
        } 

    /*
        std::cout << " * XXX_k at iter ...."<< iter << std::endl;
        for (int k = 0; k < sol.size(); ++k)
            std::cout << "  " << sol[k]  <<  std::endl;
    */

        if (convergence(sol))
        {

            var_time.pressure = var_.pressure;
            var_time.flux = var_.flux;
            var_time.L_rate = sol.tail(num_nodes_);
            return; 
        }
        
    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;
}





} //end namespace shimmer

