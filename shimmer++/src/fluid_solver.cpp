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

linearized_fluid_solver::linearized_fluid_solver(const double& tolerance, 
                        const double& dt,
                        const double& Tm,
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
                            equation_of_state *eos,
                            vector_t& sol_time)
{
    // Initialization of varibales with solution in time n; 
    vector_t press_time = sol_time.head(num_nodes_);
    vector_t flux_time  = sol_time.segment(num_nodes_, num_pipes_);

    press_ = press_time;
    flux_  = flux_time;

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

        auto mass = continuity(dt_,Tm_, press_,press_time,c2_nodes, inc_, graph_);
        auto mom  = momentum(  dt_,Tm_, flux_, flux_time, press_, press_pipes_,
                                                          c2_pipes, inc_, graph_);
        
        auto rho = eos->density_correction();
        auto bcnd = boundary(num_nodes_,num_pipes_,area_pipes, p_in, flux_,flux_ext,
                                                                rho, inlet_nodes);

        auto [LHS, rhs]= assemble(mass, mom, bcnd, graph_);


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
            sol_time.head(num_nodes_) = press_;
            sol_time.segment(num_nodes_, num_pipes_) = flux_;
            sol_time.tail(num_nodes_) = sol.tail(num_nodes_);
            return; 
        }
        
    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;
}





} //end namespace shimmer

