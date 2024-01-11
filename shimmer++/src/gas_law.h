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


namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; 
        
typedef GERG::Reducing_parameters<vector_t>   gerg_reducing_params_t;
typedef GERG::Pseudo_critical_point<vector_t> gerg_pseudo_critical_pt_t;
typedef GERG::Thermodynamic_properties_parameters gerg_thermo_params_t;
typedef GERG::Thermodynamic_properties<vector_t> gerg_thermo_props_t;


struct gerg_params
{
    gerg_reducing_params_t       reducing_params;
    gerg_pseudo_critical_pt_t    pseudo_critical_pt;
    gerg_thermo_params_t         params;  
    
    gerg_params(const gerg_reducing_params_t&  rp, 
                const gerg_pseudo_critical_pt_t& psc,
                const gerg_thermo_params_t& tp): 
                reducing_params(rp),
                pseudo_critical_pt (psc),
                params (tp) {}; 


};


gerg_thermo_props_t
equation_of_state(  const double  & temperature,
                    const vector_t& pressure,
                    const matrix_t& x,
                    const gerg_params& gerg);

} //end namespace shimmer