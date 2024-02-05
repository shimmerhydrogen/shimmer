#include <Eigen/SparseLU>
#include "../src/fluid_solver.h"
#include <iomanip>
namespace shimmer{

void
linearized_fluid_solver(const double& tolerance, 
                        const double& dt,
                        const double& Tm,
                        const incidence & inc,
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
    // old/actual/new are used for time steps 
    // previous/current/next are used for solver iterations
    size_t MAX_ITERS = 500;

    size_t num_pipes = num_edges(graph);
    size_t num_nodes = num_vertices(graph);

    // Initialization of varibales with solution in time n; 
    vector_t press_time = sol_time.head(num_nodes);
    vector_t flux_time  = sol_time.segment(num_nodes, num_pipes);

    vector_t press = press_time;
    vector_t flux  = flux_time;

	/// Underrelaxation coefficients to help convergence
    double a_G = 0.0;
    double a_p = 0.0;

    matrix_t x_nodes = build_x_nodes(graph);
    matrix_t x_pipes = inc.matrix_in().transpose() * x_nodes;

    std::cout<< "Initializing ..."<< std::endl;

    for(size_t iter = 0; iter <= MAX_ITERS; iter++)
    {    
        auto eos_nodes = equation_of_state(Tm, press, x_nodes, gerg_nodes);
        vector_t c2_nodes = eos_nodes.Z.cwiseProduct(RR_nodes) * Tm; 

        vector_t press_pipes = average(press, inc);
        auto eos_pipes = equation_of_state(Tm, press_pipes, x_pipes, gerg_pipes);
        vector_t c2_pipes = eos_pipes.Z.cwiseProduct(RR_pipes) * Tm;

        auto mass = continuity( dt,Tm, press, press_time, c2_nodes, inc, graph);
        auto mom  = momentum(dt, Tm, flux, flux_time, press, press_pipes,
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

        /*  Warning: Marco computes diff = LHS * sol - rhs, with rhs evaluated
            using the current flux (G_k). Here, the rhs is evaluated in the 
            previous step G_{k-1}. 
            res = norm(ADP*p_k -  (Rf * abs(G_k)*G_k) + Ri. * (G_k - G_n)));
   		        = norm(ADP*p_k -(2 Rf * abs(G_k)+ Ri)*G_k + (Rf *abs(G_k_1)G_k_1 + Ri * G_n))); 
        */
        // diff = lhs sol - rhs
        vector_t diff = LHS * sol - rhs; 
        auto norm_mass = diff.head(num_nodes).norm();
        auto norm_mom  = diff.segment(num_nodes, num_pipes).norm();
        auto residual  = std::max(norm_mass, norm_mom);

        vector_t press_next = a_p*press+(1.0-a_p)*sol.head(num_nodes);
        vector_t flux_next  = a_G*flux +(1.0-a_G)*sol.segment(num_nodes, num_pipes);
        
        std::cout<< "Solver at iteration k ..."<< iter << std::endl;
        if(residual < tolerance)
        {
            sol_time.head(num_nodes) = press_next;
            sol_time.segment(num_nodes, num_pipes) = flux_next;
            sol_time.tail(num_nodes) = sol.tail(num_nodes);
            
            /*
            std::cout << "sol: ("<< norm_mass<< " , " << norm_mom << ")" <<  std::endl ;
            for (int k = 0; k < sol.size(); ++k)
                std::cout << std::setprecision(16) << sol(k)<< std::endl ;
            std::cout <<  std::endl ;
            */
            return; 
        }
    
        press = press_next;
        flux  = flux_next;
    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;
}

} //end namespace shimmer