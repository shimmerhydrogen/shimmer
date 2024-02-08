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

class equation_of_state;
class linearized_fluid_solver;

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
    vector_t press_;
    vector_t press_pipes_;
    vector_t flux_;

    const incidence& inc_; 
    const infrastructure_graph& graph_;

public:

    linearized_fluid_solver(const double& tolerance, 
                        const double& dt,
                        const double& Tm,
                        const incidence & inc,
                        const infrastructure_graph& graph);

    vector_t
    boundary_velocity(equation_of_state *eos);


    bool 
    convergence(const vector_t& diff, 
                const vector_t& sol);


    void
    run(const vector_t& inlet_nodes,
        const vector_t& p_in,
        const vector_t& flux_ext,
        equation_of_state *eos,
        vector_t& sol_time);

    double temperature();         
    vector_t pressure_nodes();
    vector_t pressure_pipes();
    matrix_t x_nodes();
    matrix_t x_pipes();
    const incidence& get_incidence(); 
};





} //end namespace shimmer