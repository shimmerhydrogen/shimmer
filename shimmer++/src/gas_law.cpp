/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it  
 */

#include "../src/gas_law.h"

namespace shimmer{


gerg_thermo_props_t
equation_of_state(  const double& temperature,
                    const vector_t& pressure,
                    const matrix_t& x,
                    const gerg_params& gerg)
{
    return GERG::thermodynamic_properties(pressure, temperature, x,
                                              gerg.reducing_params,
                                              gerg.pseudo_critical_pt,
                                              gerg.params);

}

} //end namespace shimmer