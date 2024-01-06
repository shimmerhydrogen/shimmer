#include <Eigen/SparseLU>
#include "../src/fluid_solver.h"

namespace shimmer{




void
linearized_fluid_solver(const double & tolerance, const double& dt,
                        const double & Tm, const incidence& inc,
                        const infrastructure_graph& graph,
                        const vector_t& RR,
                        const gerg_params& gerg,
                        vector_t& sol_time) 
{
    size_t num_pipes = num_edges(graph);
    size_t num_nodes = num_vertices(graph);

    // Solution in time n; 
    vector_t press_old = sol_time.head(num_nodes);
    vector_t flux_old  = sol_time.segment(num_nodes, num_pipes);

    size_t MAX_ITERS = 500;

	/// Underrelaxation coefficients to help convergence
	double alfa_G = 0.0;
	double alfa_p = 0.0;

    vector_t press = press_old; 
    vector_t flux = flux_old; 

    Eigen::MatrixXd  x_nodes = build_x_nodes(graph);
    Eigen::MatrixXd  x_pipes = inc.matrix_in().transpose() * x_nodes;


    for(size_t iter = 0; iter <= MAX_ITERS; iter++){
        
        vector_t press_previous = press;
        vector_t flux_previous = flux;

        auto syst_mass = continuity(dt, Tm, press, press_old, inc, graph, x_nodes, RR, gerg);
        auto syst_mom  = momentum(dt, Tm, flux, flux_old, press, inc, graph, x_pipes, RR, gerg );
        auto [LHS, rhs]= assemble(syst_mass, syst_mom, graph);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(LHS);
        vector_t sol = solver.solve(rhs);
        vector_t press = sol.head(num_nodes);
        vector_t flux  = sol.segment(num_nodes, num_pipes);
        if(solver.info() != Eigen::Success) return;
        
        /*  Warning: Marco computes diff = lhs sol - rhs, with rhs evaluated
            using the current flux (G_k). Here, the rhs is evaluated in the 
            previous step G_{k-1}. 
            res = norm(ADP*p_k -  (Rf * abs(G_k)*G_k) + Ri. * (G_k - G_n)));
   		        = norm(ADP*p_k -(2 Rf * abs(G_k)+ Ri)*G_k + (Rf *abs(G_k_1)G_k_1 + Ri * G_n))); 
        */
        // diff = lhs sol - rhs
        vector_t diff = sol - rhs; 
        auto norm_mass = diff.head(num_nodes).norm();
        auto norm_mom  = diff.segment(num_nodes, num_pipes).norm();
        auto residual  = std::max(norm_mass, norm_mom);

        press = alfa_p * press + (1.0 - alfa_p) * press_previous;
		flux  = alfa_G * flux  + (1.0 - alfa_G) * flux_previous;

        if(residual < tolerance)
        {
            sol_time.head(num_nodes) = press;
            sol_time.segment(num_nodes, num_pipes) = flux;
            return; 
        }
    }

    std::cout << "Linearized fluid dynamics solver has NOT CONVERGED." << std::endl;

}

} //end namespace shimmer