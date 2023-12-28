#ifndef __MATLAB_GERG_functions_H
#define __MATLAB_GERG_functions_H

#include "GERG_functions.hpp"

namespace GERG
{
  // iFlag=0;
  // [Tr_b,Dr_b] = ReducingParametersGERG(Aplus'*reshape(CC_gas(:,:,1),dimn,21));
  // [Tcx_b,Dcx_b,Vcx_b]=PseudoCriticalPointGERG(Aplus'*reshape(CC_gas(:,:,1),dimn,21),dimb);
  // [Pcheck, Zm, Den]=PropertiesGERG(iFlag, pm(:,1)/1e3, Tb, Aplus'*reshape(CC_gas(:,:,1),dimn,21),dimb,Tr_b,Dr_b,Tcx_b,Dcx_b,Vcx_b);

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
  Reducing_parameters<matrix_type> reducing_parameters(const matrix_type& x)
  {
    Reducing_parameters<matrix_type> result;

    return result;
  }
}

#endif // __MATLAB_GERG_functions_H
