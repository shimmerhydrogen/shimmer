/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "solver/matrix_manipulations.h"
#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "solver/conservation_matrices.h"
#include "solver/gas_law.h"
#include "solver/variable.h"


namespace shimmer
    {

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

    vector_t c2_nodes_;
    vector_t c2_pipes_;

    vector_t T_nodes_;
    vector_t T_pipes_;

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
    impose_edge_station_model(  const vector_t& c2_nodes,
                                const vector_t& pressure_nodes,
                                const vector_t& flux,
                                sparse_matrix_t& sADP,
                                vector_t& r_scale,
                                vector_t& rhs_mom);

    pair_trip_vec_t
    momentum(const vector_t& pressure_nodes,
             const vector_t& pressure_pipes,
             const vector_t& flux,
             const vector_t& flux_old,
             const vector_t& c2_nodes,
             const vector_t& c2_pipes);


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
        equation_of_state *eos,
        size_t at_iteration = 0);

    bool check_hard_constraints(size_t step);
    void check_soft_constraints(size_t step);
    bool check_constraints(size_t step);
    bool check_hard_controls(size_t step);
    //bool check_soft_controls(size_t step);
    bool check_controls(size_t step);

    inline size_t num_nodes() const {return num_nodes_;};
    inline size_t num_pipes() const {return num_pipes_;};

    inline double temperature() const {return Tm_;};
    inline vector_t temperature_nodes() const {return T_nodes_;};
    inline vector_t temperature_pipes() const {return T_pipes_;};
    inline vector_t pressure_nodes() const {return var_.pressure;};
    inline vector_t pressure_pipes() const {return press_pipes_;};
    inline matrix_t x_nodes() const {return x_nodes_;};
    inline matrix_t x_pipes() const {return x_pipes_;};
    const incidence& get_incidence() const {return inc_;};
    inline const variable& get_variable() const{return var_;}
};





} //end namespace shimmer