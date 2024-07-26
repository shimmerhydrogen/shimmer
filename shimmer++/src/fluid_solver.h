/* This code is part of the SHIMMER project
*
* Politecnico di Torino, Dipartimento di Matematica (DISMA)
*
* Karol Cascavita (C) 2023
* karol.cascavita@polito.it
*/

#pragma once

#include "../src/matrix_manipulations.h"
#include "../src/infrastructure_graph.h"
#include "../src/incidence_matrix.h"
#include "../src/conservation_matrices.h"
#include "../src/gas_law.h"
#include "../src/temporal.h"

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>;
using triplet_t = Eigen::Triplet<double>;
using pair_trip_vec_t = std::pair<std::vector<triplet_t>, vector_t>;

class equation_of_state;
//class linearized_fluid_solver;

class linearized_fluid_solver
{
    bool is_unsteady_;
    size_t MAX_ITERS_;
    size_t num_pipes_;
    size_t num_nodes_;
    size_t at_step_;

    double tolerance_;
    double a_G_;
    double a_p_;
    double dt_;
    double Tm_;

    matrix_t x_nodes_;
    matrix_t x_pipes_;
    vector_t press_pipes_;
    vector_t mu_;

    variable var_;
    const incidence& inc_;
    const infrastructure_graph& graph_;

public:

    linearized_fluid_solver(size_t at_step, const bool& is_unsteady,
                        double tolerance,
                        double dt,
                        double Tm,
                        const vector_t& mu,
                        const incidence & inc,
                        const infrastructure_graph& graph);


    pair_trip_vec_t
    continuity(const vector_t& pressure_old,
            const vector_t& c2);

    void
    control_stations(std::vector<triplet_t>& triplets_mom,
                     vector_t& rhs_mom);

    pair_trip_vec_t
    momentum(
             const vector_t& pressure_nodes,
             const vector_t& pressure_pipes,
             const vector_t& flux,
             const vector_t& flux_old,
             const vector_t& c2);


    pair_trip_vec_t
    boundary(const vector_t& area_pipes,
            const vector_t& flux,
            equation_of_state *eos);


    std::pair<sparse_matrix_t, vector_t>
    assemble(const pair_trip_vec_t& lhs_rhs_mass,
             const pair_trip_vec_t& lhs_rhs_mom,
             const pair_trip_vec_t& lhs_rhs_bnd);


    bool
    convergence(const vector_t& sol);


    bool
    run(const vector_t& area_pipes,
        //const vector_t& flux_ext,
        const variable& var_guess,
        const variable& var_time,
        equation_of_state *eos
        );
    bool check_hard_constraints(size_t step);
    void check_soft_constraints(size_t step);
    bool check_constraints(size_t step);
    bool check_edge_constraints(size_t step);

    double temperature();
    vector_t pressure_nodes();
    vector_t pressure_pipes();
    matrix_t x_nodes();
    matrix_t x_pipes();
    const incidence& get_incidence();
    inline const variable& get_variable(){return var_;}
};





} //end namespace shimmer