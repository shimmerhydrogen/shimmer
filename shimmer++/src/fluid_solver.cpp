#include <Eigen/SparseLU>
#include "../src/fluid_solver.h"
#include <iomanip>
namespace shimmer{


linearized_fluid_solver::linearized_fluid_solver(const double& tolerance, 
                    const double& dt,
                    const double& Tm_,
                    const incidence & inc,
                    const infrastructure_graph& graph): 
                    tolerance_(tolerance), dt_(dt), Tm_(Tm_)  
{
    MAX_ITERS_ = 500;

    num_pipes_ = num_edges(graph);
    num_nodes_ = num_vertices(graph);

    /// Underrelaxation coefficients to help convergence
    a_G_ = 0.0;
    a_p_ = 0.0;

    x_nodes_ = build_x_nodes(graph);
    x_pipes_ = inc.matrix_in().transpose() * x_nodes_;
}


bool 
linearized_fluid_solver::convergence(const vector_t& diff, 
                const vector_t& sol,
                vector_t& flux,
                vector_t& press) 
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

    vector_t press_next = a_p_*press+(1.0-a_p_)*sol.head(num_nodes_);
    vector_t flux_next  = a_G_*flux +(1.0-a_G_)*sol.segment(num_nodes_, num_pipes_);

    press = press_next;
    flux  = flux_next;
    
    if(residual < tolerance_)
    {
        /*
        std::cout << "sol: ("<< norm_mass<< " , " << norm_mom << ")" <<  std::endl ;
        for (int k = 0; k < sol.size(); ++k)
            std::cout << std::setprecision(16) << sol(k)<< std::endl ;
        std::cout <<  std::endl ;
        */
        return true; 
    }


    return false;
}

            
void
linearized_fluid_solver::compute(const incidence & inc,
        const infrastructure_graph& graph,
        const vector_t& inlet_nodes,
        const double& p_in,
        const vector_t& flux_ext,
        const vector_t& RR_nodes,
        const vector_t& RR_pipes,
        const vector_t& molar_mass,
        const gerg_params& gerg_nodes,
        const gerg_params& gerg_pipes,
        vector_t& sol_time)
{
    // Initialization of varibales with solution in time n; 
    vector_t press_time = sol_time.head(num_nodes_);
    vector_t flux_time  = sol_time.segment(num_nodes_, num_pipes_);

    vector_t press = press_time;
    vector_t flux  = flux_time;

    std::cout<< "Initializing ..."<< std::endl;

    for(size_t iter = 0; iter <= MAX_ITERS_; iter++)
    {   

        std::cout<< "Solver at iteration k ..."<< iter << std::endl;
    
        //auto [c2_nodes, c2_pipes] = speed_of_sound(Tm_, press, press_pipes, this); 
        auto eos_nodes = equation_of_state(Tm_, press, x_nodes_, gerg_nodes);
        vector_t c2_nodes = eos_nodes.Z.cwiseProduct(RR_nodes) * Tm_; 

        vector_t press_pipes = average(press, inc);
        auto eos_pipes = equation_of_state(Tm_, press_pipes, x_pipes_, gerg_pipes);
        vector_t c2_pipes = eos_pipes.Z.cwiseProduct(RR_pipes) * Tm_;

        auto mass = continuity(dt_,Tm_, press, press_time, c2_nodes, inc, graph);
        auto mom  = momentum(dt_, Tm_, flux, flux_time, press, press_pipes,
                                                    c2_pipes, inc, graph);
        auto bcnd = boundary(p_in, flux, flux_ext, molar_mass, 
                                    inc, graph, inlet_nodes, eos_pipes);

        auto [LHS, rhs]= assemble(mass, mom, bcnd, graph);


        Eigen::SparseLU<sparse_matrix_t> solver;
        solver.compute(LHS);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error factorizing LHS" <<std::endl;
            exit(1);
        }
        vector_t sol = solver.solve(rhs);
        if(solver.info() != Eigen::Success) {
            std::cout << "Error solving system" <<std::endl;
            exit(1);
        } 

        vector_t diff = LHS * sol - rhs; 

        if (convergence(diff, sol, flux, press))
        {
            sol_time.head(num_nodes_) = press;
            sol_time.segment(num_nodes_, num_pipes_) = flux;
            sol_time.tail(num_nodes_) = sol.tail(num_nodes_);
            return; 
        }
        
    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;
}


} //end namespace shimmer

