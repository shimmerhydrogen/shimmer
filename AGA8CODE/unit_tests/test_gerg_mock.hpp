#ifndef __test_gerg_mock_H
#define __test_gerg_mock_H

#include "Eigen/Eigen"
#include "shimmer_gerg_data.hpp"

namespace shimmer_gerg
{
  namespace test
  {
    namespace gerg_mock
    {
      // *********************************************************
      double tolerance() { return 1.0e-12; }
      // *********************************************************
      auto x()
      {
        std::vector<double> x(21, 0.0);
        x[0] = 0.8;
        x[2] = 0.2;

        return x;
      }
      // *********************************************************
      auto x_mat()
      {
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(1,21);
        x(0,0) = 0.8;
        x(0,2) = 0.2;

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
        std::vector<double> mol_frac(2);
        mol_frac[0] = 0.8;
        mol_frac[1] = 0.2;

        return mol_frac;
      }
      // *********************************************************
      auto comps()
      {
        std::vector<std::string> comps = { "Methane", "Carbon dioxide" };

        return comps;
      }
      // *********************************************************
      auto reducing_parameters()
      {
        gerg_data::Reducing_parameters<double> reducing_parameters;

        reducing_parameters.Tr = 2.082743529542123e+02;
        reducing_parameters.Dr = 1.022348540995680e+01; // mol/L * 1e3 1/m^3 = mol/m^3

        return reducing_parameters;
      }
      // *********************************************************
      auto pseudo_critical_point()
      {
        gerg_data::Pseudo_critical_point<double> pseudo_critical_point;

        pseudo_critical_point.Tcx = 2.132768400000000e+02;
        pseudo_critical_point.Vcx = 9.772414504557397e-02;
        pseudo_critical_point.Dcx = 1.023288563469802e+01;

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
      auto alpha_r()
      {
        return -1.146681132712004e-01;
      }
      // *********************************************************
      auto alpha_0()
      {
         return 2.539691191428680e+00;
      }
      // *********************************************************
      auto delta()
      {
        return 2.312348728992204e-01;
      }
      // *********************************************************
      auto tau()
      {
        return 7.104702471574701e-01;
      }
      // *********************************************************
      auto thermodynamic_properties_parameters()
      {
        gerg_data::Thermodynamic_properties_parameters<double> parameters;
        parameters.Type = gerg_data::Thermodynamic_properties_parameters<double>::Types::Gas_phase;
        parameters.T = T();
        parameters.P = P();

        return parameters;
      }
      // *********************************************************
      auto thermodynamic_properties()
      {
        gerg_data::Thermodynamic_properties<double> thermodynamic_properties;
        thermodynamic_properties.D = 4.853135975758211e+01;
        thermodynamic_properties.Z = 4.312568035613069e+01;
        thermodynamic_properties.P = 5.101324999999521e+06;
        thermodynamic_properties.gamma = 1.1714382285678768e+00;
        thermodynamic_properties.Cp = 5.2298638227860138e+01;
        thermodynamic_properties.Cv = 4.4644810927672225e+01;
        thermodynamic_properties.dPdD = 5.3312067323673121e+05;
        thermodynamic_properties.dPdT = 5.7257133461214125e+03;

        return thermodynamic_properties;
      }
      // *********************************************************
      auto thermodynamic_properties_mat()
      {
        gerg_data::Thermodynamic_properties<Eigen::VectorXd> thermodynamic_properties;
        thermodynamic_properties.D = Eigen::VectorXd::Zero(1);
        thermodynamic_properties.Z = Eigen::VectorXd::Zero(1);
        thermodynamic_properties.P = Eigen::VectorXd::Zero(1);
        thermodynamic_properties.gamma= Eigen::VectorXd::Zero(1);
        thermodynamic_properties.Cp   = Eigen::VectorXd::Zero(1);
        thermodynamic_properties.Cv   = Eigen::VectorXd::Zero(1);
        thermodynamic_properties.dPdD = Eigen::VectorXd::Zero(1);
        thermodynamic_properties.dPdT = Eigen::VectorXd::Zero(1);

        thermodynamic_properties.D[0] = 4.853135975758211e+01;
        thermodynamic_properties.Z[0] = 4.312568035613069e+01;
        thermodynamic_properties.P[0] = 5.101324999999521e+06;
        thermodynamic_properties.gamma[0]= 1.1714382285678768e+00;
        thermodynamic_properties.Cp[0]   = 5.2298638227860138e+01;
        thermodynamic_properties.Cv[0]   = 4.4644810927672225e+01;
        thermodynamic_properties.dPdD[0] = 5.3312067323673121e+05;
        thermodynamic_properties.dPdT[0] = 5.7257133461214125e+03;

        return thermodynamic_properties;
      }
      // *********************************************************
    }
  }
}

#endif // __test_temp_H
