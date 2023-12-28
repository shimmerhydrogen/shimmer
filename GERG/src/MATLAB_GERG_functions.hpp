#ifndef __MATLAB_GERG_functions_H
#define __MATLAB_GERG_functions_H

#include "Matlab_interface.hpp"

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
    const Matlab_interface& matlab = Matlab_interface::get_instance();

    matlab::data::ArrayFactory factory;

    const matlab::data::TypedArray<double> x_to_matlab = matlab.matrix_to_matlab(factory,
                                                                                 x);

    std::cout<< "Test\n";
    for (unsigned int r = 0; r < x.rows(); r++)
    {
      for (unsigned int c = 0; c < x.cols(); c++)
        std::cout<< x_to_matlab[r][c]<< ", ";
      std::cout<< std::endl;
    }

    Matlab_interface::get_instance().reducing_parameters(x_to_matlab);

    Reducing_parameters<matrix_type> result;

    return result;
  }
}

#endif // __MATLAB_GERG_functions_H
