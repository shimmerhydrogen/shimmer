#ifndef __test_temp_2_H
#define __test_temp_2_H

#include <iostream>
#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/cpp/deriv_adapter.hpp"
#include "teqp/cpp/derivs.hpp"
#include "teqp/models/GERG/GERG.hpp"

#include "shimmer_teqp_utilities.hpp"

#include "test_utilities.hpp"

namespace shimmer_teqp
{
  namespace test
  {
    // *********************************************************
    int test_temp_2(int , char **)
    {
      using namespace teqp::cppinterface;
      using namespace teqp::cppinterface::adapter;

      constexpr double tolerance = std::numeric_limits<double>::epsilon();
      Eigen::ArrayXd x(21);
      x<< 0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

      const auto filter = shimmer_teqp::utilities::filter_components(x,
                                                                     tolerance);

      const auto& names = teqp::GERG2008::component_names;
      const auto comps = shimmer_teqp::utilities::filter_collection(filter,
                                                                    names);
      const auto mol_frac = shimmer_teqp::utilities::filter_collection(filter,
                                                                       x);

      // std::cout<< "Comps: "<< comps<< std::endl;
      // std::cout<< "Mol_frac: "<< mol_frac<< std::endl;

      auto model = teqp::GERG2008::GERG2008ResidualModel(comps);

      const double T = 300;
      const double rho_molar = 3;

      const auto R = model.R(mol_frac);
      std::cout<< "R: " << R << std::endl;
      const auto Tc = model.red.get_Tcvec();
      const auto Vc = model.red.get_vcvec();
      const auto Tr = model.red.get_Tr(mol_frac);
      const auto rho_r = model.red.get_rhor(mol_frac);


      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_2_H
