#ifndef __test_gerg_mock_H
#define __test_gerg_mock_H

#include "Eigen/Eigen"
#include "shimmer_gerg_data.hpp"

namespace shimmer_teqp
{
  namespace test
  {
    namespace gerg_mock
    {
      // *********************************************************
      double tolerance() { return 1.0e-14; }
      // *********************************************************
      auto x()
      {
        Eigen::ArrayXd x(21);
        x<< 0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

        return x;
      }
      // *********************************************************
      auto filter()
      {
        std::vector<unsigned int> filter = { 0, 2 };

        return filter;
      }
      // *********************************************************
      auto mol_frac()
      {
        Eigen::ArrayXd mol_frac(2);
        mol_frac<< 0.8, 0.2;

        return mol_frac;
      }
      // *********************************************************
      auto comps()
      {
        std::vector<std::string> comps = { "methane", "carbondioxide" };

        return comps;
      }
      // *********************************************************
      auto reducing_parameters()
      {
        gerg_data::Reducing_parameters<double> reducing_parameters;

        reducing_parameters.Tr = 2.082743529542123e+02;
        reducing_parameters.Dr = 1.022348540995680e+01 * 1.0e+3;

        return reducing_parameters;
      }
      // *********************************************************
      auto pseudo_critical_point()
      {
        gerg_data::Pseudo_critical_point<double> pseudo_critical_point;

        pseudo_critical_point.Tcx = 2.132768400000000e+02;
        pseudo_critical_point.Vcx = 9.772414504557397e-02 * 1.0e-3;
        pseudo_critical_point.Dcx = 1.023288563469802e+01 * 1.0e+3;

        return pseudo_critical_point;
      }
      // *********************************************************
      auto T()
      {
        return 2.931500000000000e+02;
      }
      // *********************************************************
      auto P()
      {
        return 5.101325000000000e+06;
      }
      // *********************************************************
      auto thermodynamic_properties_parameters()
      {
        gerg_data::Thermodynamic_properties_parameters parameters;
        parameters.Type = gerg_data::Thermodynamic_properties_parameters::Types::Gas_phase;

        return parameters;
      }
      // *********************************************************
      auto thermodynamic_properties()
      {
        gerg_data::Thermodynamic_properties<double> thermodynamic_properties;
        thermodynamic_properties.D = 4.853135975758211e+01;
        thermodynamic_properties.Z = 4.312568035613069e+01;
        thermodynamic_properties.P1 = 5.101324999999521e+06;
        thermodynamic_properties.gamma = 3.571899838350778e+00;

        return thermodynamic_properties;
      }
      // *********************************************************
    }
  }
}

#endif // __test_temp_H
