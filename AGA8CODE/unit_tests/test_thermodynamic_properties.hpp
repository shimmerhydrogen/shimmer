#ifndef __test_thermodynamic_properties_H
#define __test_thermodynamic_properties_H

#include "Eigen/Eigen"
#include "shimmer_gerg_functions.hpp"
#include "shimmer_teqp_utilities.hpp"

#include "test_gerg_mock.hpp"
#include "test_utilities.hpp"

namespace shimmer_teqp
{
  namespace test
  {
    // *********************************************************
    int test_thermodynamic_properties(int , char **)
    {
      using namespace shimmer_teqp::utilities;

      const auto tolerance = gerg_mock::tolerance();
      const auto comps = gerg_mock::comps();
      const auto mol_frac = gerg_mock::mol_frac();

      const auto x = gerg_mock::x();

      shimmer_teqp::gerg_functions::setup_GERG();

      const auto input_properties = gerg_mock::thermodynamic_properties_parameters();
      const auto thermodynamic_properties = shimmer_teqp::gerg_functions::thermodynamic_properties(x,
                                                                                                   input_properties,
                                                                                                   tolerance);
      const auto expected_thermodynamic_properties = gerg_mock::thermodynamic_properties();

      std::cout.precision(16);
      std::cout<< "ASK to MARCO the differences on MATLAB computation, as parameter alpha_r(:,2, 2) is 0.0 in MATLAB"<< std::endl;
      std::cout<< std::scientific<< "\t C++ gamma "<< thermodynamic_properties.gamma<< std::endl;
      std::cout<< std::scientific<< "\t C++ dPdD "<< thermodynamic_properties.dPdD<< std::endl;
      std::cout<< std::scientific<< "\t C++ dPdT "<< thermodynamic_properties.dPdT<< std::endl;
      std::cout<< std::scientific<< "\t C++ Cp "<< thermodynamic_properties.Cp<< std::endl;
      std::cout<< std::scientific<< "\t C++ Cv "<< thermodynamic_properties.Cv<< std::endl;

      std::cout<< std::scientific<< "\t MATLAB gamma "<< expected_thermodynamic_properties.gamma<< std::endl;
      std::cout<< std::scientific<< "\t MATLAB dPdD "<< expected_thermodynamic_properties.dPdD<< std::endl;
      std::cout<< std::scientific<< "\t MATLAB dPdT "<< expected_thermodynamic_properties.dPdT<< std::endl;
      std::cout<< std::scientific<< "\t MATLAB Cp "<< expected_thermodynamic_properties.Cp<< std::endl;
      std::cout<< std::scientific<< "\t MATLAB Cv "<< expected_thermodynamic_properties.Cv<< std::endl;

      ASSERT_DOUBLE_EQ_TOL(expected_thermodynamic_properties.D,
                           thermodynamic_properties.D,
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(expected_thermodynamic_properties.Z,
                           thermodynamic_properties.Z,
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(expected_thermodynamic_properties.P,
                           thermodynamic_properties.P,
                           tolerance);
      // ASSERT_DOUBLE_EQ_TOL(expected_thermodynamic_properties.gamma,
      //                      thermodynamic_properties.gamma,
      //                      tolerance);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
