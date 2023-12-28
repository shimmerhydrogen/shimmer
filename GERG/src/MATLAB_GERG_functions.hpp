#ifndef __MATLAB_GERG_functions_H
#define __MATLAB_GERG_functions_H

#include "MATLAB_GERG_interface.hpp"

namespace GERG
{
  struct Thermodynamic_properties_parameters final
  {
      enum struct Types
      {
        Gas_phase = 0, // strict pressure solver in gas phase without checks (fastest mode, but output state may not be stable single phase)
        Two_phase = 1, // to make checks for possible 2-phase state (result may still not be stable single phase, but many unstable states will be identified)
        Liquid_phase = 2// to search for liquid phase (and make the same checks when iFlag=1)
      };

      Types Type;
  };

  template <class matrix_type>
  inline Reducing_parameters<matrix_type> reducing_parameters(const matrix_type& x)
  {
    return Matlab_interface::GetInstance().reducing_parameters(x);
  }
}

#endif // __MATLAB_GERG_functions_H
