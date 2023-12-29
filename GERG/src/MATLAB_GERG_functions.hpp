#ifndef __MATLAB_GERG_functions_H
#define __MATLAB_GERG_functions_H

#include "GERG_functions.hpp"
#include "Matlab_interface.hpp"

namespace GERG
{
  // *********************************************************
  struct Thermodynamic_properties_parameters final
  {
      enum struct Types
      {
        Gas_phase = 0, ///< strict pressure solver in gas phase without checks (fastest mode, but output state may not be stable single phase)
        Two_phase = 1, ///< to make checks for possible 2-phase state (result may still not be stable single phase, but many unstable states will be identified)
        Liquid_phase = 2 ///< to search for liquid phase (and make the same checks when iFlag=1)
      };

      Types Type;
  };
  // *********************************************************
  template <class matrix_type>
  inline Reducing_parameters<matrix_type> reducing_parameters(const matrix_type& x)
  {
    const Matlab_interface& matlab = Matlab_interface::get_instance();

    matlab::data::ArrayFactory factory;

    const matlab::data::TypedArray<double> x_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                            x);

    const std::vector<matlab::data::Array> reducing_parameters = matlab.reducing_parameters(x_to_matlab);


    Reducing_parameters<matrix_type> result;
    result.Tr = Matlab_interface::matlab_to_matrix<matrix_type>(reducing_parameters.at(0));
    result.Dr = Matlab_interface::matlab_to_matrix<matrix_type>(reducing_parameters.at(1));

    return result;
  }
  // *********************************************************
  template <class matrix_type>
  Pseudo_critical_point<matrix_type> pseudo_critical_point(const matrix_type& x,
                                                           const unsigned int dimn)
  {
    const Matlab_interface& matlab = Matlab_interface::get_instance();

    matlab::data::ArrayFactory factory;

    const matlab::data::TypedArray<double> x_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                            x);
    const matlab::data::TypedArray<double> dimn_to_matlab = factory.createScalar(static_cast<double>(dimn));

    const std::vector<matlab::data::Array> pseudo_critical_point = matlab.pseudo_critical_point(x_to_matlab,
                                                                                                dimn_to_matlab);

    Pseudo_critical_point<matrix_type> result;
    result.Tcx = Matlab_interface::matlab_to_matrix<matrix_type>(pseudo_critical_point.at(0));
    result.Dcx = Matlab_interface::matlab_to_matrix<matrix_type>(pseudo_critical_point.at(1));
    result.Vcx = Matlab_interface::matlab_to_matrix<matrix_type>(pseudo_critical_point.at(2));

    return result;
  }
  // *********************************************************
  template <class matrix_type>
  Thermodynamic_properties<matrix_type> thermodynamic_properties(const matrix_type& P,
                                                                 const matrix_type& T,
                                                                 const matrix_type& x,
                                                                 const unsigned int dimn,
                                                                 const matrix_type& Tr,
                                                                 const matrix_type& Dr,
                                                                 const matrix_type& Tcx,
                                                                 const matrix_type& Dcx,
                                                                 const matrix_type& Vcx,
                                                                 const Thermodynamic_properties_parameters& parameters)
  {
    const Matlab_interface& matlab = Matlab_interface::get_instance();

    matlab::data::ArrayFactory factory;

    const matlab::data::TypedArray<double> P_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                            P);
    const matlab::data::TypedArray<double> T_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                            T);
    const matlab::data::TypedArray<double> x_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                            x);
    const matlab::data::TypedArray<double> dimn_to_matlab = factory.createScalar(static_cast<double>(dimn));
    const matlab::data::TypedArray<double> Tr_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                             Tr);
    const matlab::data::TypedArray<double> Dr_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                             Dr);
    const matlab::data::TypedArray<double> Tcx_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                              Tcx);
    const matlab::data::TypedArray<double> Dcx_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                              Dcx);
    const matlab::data::TypedArray<double> Vcx_to_matlab = Matlab_interface::matrix_to_matlab(factory,
                                                                                              Vcx);
    const matlab::data::TypedArray<double> iFlag_to_matlab = factory.createScalar(static_cast<double>(parameters.Type));

    const std::vector<matlab::data::Array> thermodynamic_properties = matlab.thermodynamic_properties(P_to_matlab,
                                                                                                      T_to_matlab,
                                                                                                      x_to_matlab,
                                                                                                      dimn_to_matlab,
                                                                                                      Tr_to_matlab,
                                                                                                      Dr_to_matlab,
                                                                                                      Tcx_to_matlab,
                                                                                                      Dcx_to_matlab,
                                                                                                      Vcx_to_matlab,
                                                                                                      iFlag_to_matlab);

    Thermodynamic_properties<matrix_type> result;
    result.P1 = Matlab_interface::matlab_to_matrix<matrix_type>(thermodynamic_properties.at(0));
    result.Z = Matlab_interface::matlab_to_matrix<matrix_type>(thermodynamic_properties.at(1));
    result.D = Matlab_interface::matlab_to_matrix<matrix_type>(thermodynamic_properties.at(2));
    result.gamma = Matlab_interface::matlab_to_matrix<matrix_type>(thermodynamic_properties.at(3));

    return result;
  }
  // *********************************************************
}

#endif // __MATLAB_GERG_functions_H
