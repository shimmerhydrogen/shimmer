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

#include <Eigen/Dense>
#ifdef HAVE_MATLAB_GERG
#include "MATLAB_GERG_functions.hpp"
#endif /* HAVE_MATLAB_GERG */
#include "solver/incidence_matrix.h"
#include "infrastructure_graph.h"
#include "solver/fluid_solver.h"
#include "gerg/shimmer_gerg_functions.hpp"
#include "gerg/shimmer_gerg_utilities.hpp"

namespace shimmer{

class equation_of_state;
class linearized_fluid_solver;
class gerg;

#define NUM_GASES 21

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; 

matrix_t build_x_nodes(const infrastructure_graph& g);

class equation_of_state
{

protected:
    vector_t mm_nodes_;
    vector_t mm_pipes_;

    double Runiversal_ = 8314.4598;  
    vector_t Rspecific_;

public:
    equation_of_state(){};
    vector_t density(linearized_fluid_solver *) const;
    vector_t density_nodes(linearized_fluid_solver *) const;
    vector_t compute_R(const vector_t& molar_mass) const;
    void mixture_molar_mass(const infrastructure_graph&, const incidence&);

    virtual void initialization(linearized_fluid_solver *) = 0; 
    virtual void mixture_molar_mass(const matrix_t&, const matrix_t&) = 0;
    virtual std::pair<matrix_t, matrix_t> molarfrac_2_massfrac(const infrastructure_graph&, const incidence&) = 0;
    virtual std::pair<matrix_t, matrix_t> massfrac_2_molarfrac(const matrix_t&, const matrix_t&) = 0;

    virtual std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *) const = 0;
    inline vector_t mm_nodes() const {return mm_nodes_;};
    inline vector_t mm_pipes() const {return mm_pipes_;};
};


class papay: public equation_of_state
{
    double T_cr_;
    double p_cr_;

public:
    papay();
    
    using equation_of_state::mixture_molar_mass;
    void initialization(linearized_fluid_solver *lfs); 
    void mixture_molar_mass(const matrix_t&, const matrix_t&);
    std::pair<matrix_t, matrix_t>
    molarfrac_2_massfrac(const infrastructure_graph&, const incidence&);
    std::pair<matrix_t, matrix_t> 
    massfrac_2_molarfrac(const matrix_t&, const matrix_t&);
    vector_t
    compute_Z(double temperature, const vector_t& pressure) const;

    std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *lfs) const;
};


#ifdef HAVE_MATLAB_GERG 
typedef GERG::Reducing_parameters<matrix_t>   gerg_reducing_params_t;
typedef GERG::Pseudo_critical_point<matrix_t> gerg_pseudo_critical_pt_t;
typedef GERG::Thermodynamic_properties_parameters gerg_thermo_params_t;
typedef GERG::Thermodynamic_properties<vector_t> gerg_thermo_props_t;

struct gerg_params
{
    gerg_reducing_params_t       reducing_params;
    gerg_pseudo_critical_pt_t    pseudo_critical_pt;
    gerg_thermo_params_t         params;  
    
    gerg_params(const gerg_reducing_params_t&  rp, 
                const gerg_pseudo_critical_pt_t& psc,
                const gerg_thermo_params_t& tp); 

    gerg_params();
};

class gerg: public equation_of_state
{
    gerg_params params_nodes_;
    gerg_params params_pipes_;

    vector_t mm_component;
    std::vector<int> gas_name_;
    gerg_params compute_params(const matrix_t& x);

public:

    gerg();

    using equation_of_state::mixture_molar_mass;
    void initialization(linearized_fluid_solver *lfs); 
    void mixture_molar_mass(const matrix_t& x_nodes, const matrix_t& x_pipes);
    std::pair<matrix_t, matrix_t>
    molarfrac_2_massfrac(const infrastructure_graph&, const incidence&);
    std::pair<matrix_t, matrix_t> 
    massfrac_2_molarfrac(const matrix_t&, const matrix_t&);

    vector_t
    compute_Z(const double  & temperature,
              const vector_t& pressure,
              const matrix_t& x,
              const gerg_params& gp) const;


    std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *lfs) const;
};
#endif /* HAVE_MATLAB_GERG  */


using namespace shimmer_gerg;
typedef gerg_data::Thermodynamic_properties_parameters<double> gerg_aga_thermo_params_t;
typedef gerg_data::Thermodynamic_properties<vector_t> gerg_aga_thermo_props_t;

class gerg_aga: public equation_of_state
{
    double tolerance_;
    vector_t mm_component;

public:

    gerg_aga();

    using equation_of_state::mixture_molar_mass;
    void initialization(linearized_fluid_solver *lfs); 
    void mixture_molar_mass(const matrix_t& x_nodes, const matrix_t& x_pipes);
    std::pair<matrix_t, matrix_t> 
    molarfrac_2_massfrac(const infrastructure_graph&, const incidence&);
    std::pair<matrix_t, matrix_t> 
    massfrac_2_molarfrac(const matrix_t&, const matrix_t&);

    vector_t
    compute_Z(const vector_t& temperature,
            const vector_t& pressure,
            const matrix_t& x) const;


    std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *lfs) const;
};


} //end namespace shimmer
