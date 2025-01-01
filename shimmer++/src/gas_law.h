/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it  
 */

#pragma once

#include <Eigen/Dense>
#include "MATLAB_GERG_functions.hpp"
#include "../src/fluid_solver.h"

namespace shimmer{

class equation_of_state;
class linearized_fluid_solver;
class gerg;



using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; 


class equation_of_state
{

protected:
    vector_t density_;
    vector_t R_nodes_;
    vector_t R_pipes_;

    vector_t mm_nodes_;
    vector_t mm_pipes_;

    double Runiversal_ = 8314.4598;  
    vector_t Rspecific_;

public:
    equation_of_state(){};
    void compute_density(linearized_fluid_solver *, const vector_t&);

    virtual void initialization(linearized_fluid_solver *) = 0; 
    virtual void compute_molar_mass(const matrix_t&, const matrix_t&) = 0;
    virtual std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *) = 0;

    const vector_t& density(); 
};


class papay: public equation_of_state
{
    double T_cr_;
    double p_cr_;

public:
    papay();
    
    void initialization(linearized_fluid_solver *lfs); 

    void compute_molar_mass(const matrix_t& y_nodes, const matrix_t& y_pipes);

    vector_t
    compute(double temperature, const vector_t& pressure);

    std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *lfs);
};


        
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

    vector_t mmi_gerg;
    gerg_params compute_params(const matrix_t& x);

public:

    gerg();

    void compute_molar_mass(const matrix_t& y_nodes, const matrix_t& y_pipes);

    void initialization(linearized_fluid_solver *lfs); 

    vector_t
    compute(const double  & temperature,
            const vector_t& pressure,
            const matrix_t& x,
            const gerg_params& gp);


    std::pair<vector_t, vector_t>
    speed_of_sound(linearized_fluid_solver *lfs);
};





} //end namespace shimmer