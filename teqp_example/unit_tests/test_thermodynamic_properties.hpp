#ifndef __test_thermodynamic_properties_H
#define __test_thermodynamic_properties_H

#include "Eigen/Eigen"
#include "shimmer_gerg_functions.hpp"
#include "shimmer_teqp_utilities.hpp"

#include "teqp/models/GERG/GERG.hpp"
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

      auto model = teqp::GERG2008::GERG2008ResidualModel(comps);
      const auto R = model.R(mol_frac);

      const auto T = gerg_mock::T();
      const auto D = gerg_mock::thermodynamic_properties().D;
      const auto alpha_r = model.alphar(T,
                                        D,
                                        mol_frac);
      const auto Z = 1 + alpha_r;
      const auto P = D * R * T * Z;

      std::cout.precision(4);
      std::cout<< std::scientific<< "D "<< D<< std::endl;
      std::cout<< std::scientific<< "T "<< T<< std::endl;
      std::cout<< std::scientific<< "P "<< P<< std::endl;
      std::cout<< std::scientific<< "Z "<< Z<< std::endl;
      std::cout<< std::scientific<< "a "<< alpha_r<< std::endl;

      ASSERT_DOUBLE_EQ_TOL(gerg_mock::alpha_r(),
                           alpha_r,
                           tolerance);

      const auto expected_thermodynamic_properties = gerg_mock::thermodynamic_properties();

      const auto thermodynamic_properties = gerg_functions::thermodynamic_properties(mol_frac,
                                                                                     gerg_mock::tolerance());

      ASSERT_DOUBLE_EQ_TOL(expected_thermodynamic_properties.D,
                           thermodynamic_properties.D,
                           tolerance);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
