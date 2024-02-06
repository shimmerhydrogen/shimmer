/* This code is part of the SHIMMER project
*
* Politecnico di Torino, Dipartimento di Matematica (DISMA)
* 
* Karol Cascavita (C) 2023
* karol.cascavita@polito.it  
*/

#pragma once

#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/conservation_matrices.h"
#include "../src/assemble.h"
#include "../src/gas_law.h"

namespace shimmer{


class linearized_fluid_solver
{
    size_t MAX_ITERS_;
    size_t num_pipes_;
    size_t num_nodes_;

    double tolerance_;
    double a_G_;
    double a_p_;
    double dt_; 
    double Tm_;

    matrix_t x_nodes_;
    matrix_t x_pipes_;

public:

    linearized_fluid_solver(const double& tolerance, 
                        const double& dt,
                        const double& Tm_,
                        const incidence & inc,
                        const infrastructure_graph& graph); 


    bool convergence(const vector_t& diff, 
                     const vector_t& sol,
                     vector_t& flux,
                     vector_t& press); 


    void
    compute(const incidence & inc,
            const infrastructure_graph& graph,
            const vector_t& inlet_nodes,
            const double& p_in,
            const vector_t& flux_ext,
            const vector_t& RR_nodes,
            const vector_t& RR_pipes,
            const vector_t& molar_mass,
            const gerg_params& gerg_nodes,
            const gerg_params& gerg_pipes,
            vector_t& sol_time);
};

} //end namespace shimmer