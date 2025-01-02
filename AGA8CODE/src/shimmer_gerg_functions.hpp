#ifndef __SHIMMER_TEQP_FUNCTIONS_HPP__
#define __SHIMMER_TEQP_FUNCTIONS_HPP__

#include "GERG2008.h"
#include "shimmer_gerg_data.hpp"
#include <stdexcept>

namespace shimmer_teqp
{
  namespace gerg_functions
  {
    static const unsigned int GERG_num_componets = 21;
    // *********************************************************
    template<typename x_type,
             typename value_type>
    auto create_x_GERG(const x_type& x,
                       const value_type& tolerance)
    {
      using size_type = typename x_type::size_type;

      assert(static_cast<size_type>(x.size()) == GERG_num_componets);

      std::vector<double> x_GERG(GERG_num_componets + 1);

      x_GERG[0] = 0.0;
      for (unsigned int i = 0; i < GERG_num_componets; ++i)
        x_GERG[i + 1] = abs(x[i]) > tolerance ? x[i] : 0.0;

      return x_GERG;
    }
    // *********************************************************
    inline void setup_GERG() noexcept
    {
      SetupGERG();
    }
    // *********************************************************
    template<typename x_type,
             typename value_type>
    auto reducing_parameters(const x_type& x,
                             const value_type& tolerance)
    {
      using size_type = typename x_type::size_type;

      assert(static_cast<size_type>(x.size()) == GERG_num_componets);

      const auto x_GERG = create_x_GERG(x, tolerance);
      gerg_data::Reducing_parameters<value_type> reducing_parameters;

      reducing_parameters.Tr = value_type();
      reducing_parameters.Dr = value_type();

      ReducingParametersGERG(x_GERG,
                             reducing_parameters.Tr,
                             reducing_parameters.Dr);

      return reducing_parameters;
    }
    // *********************************************************
    template<typename x_type,
             typename value_type>
    auto pseudo_critical_point(const x_type& x,
                               const value_type& tolerance)
    {
      using size_type = typename x_type::size_type;

      assert(static_cast<size_type>(x.size()) == GERG_num_componets);

      gerg_data::Pseudo_critical_point<value_type> pseudo_critical_point;

      pseudo_critical_point.Tcx = value_type();
      pseudo_critical_point.Vcx = value_type();
      pseudo_critical_point.Dcx = value_type();

      const auto x_GERG = create_x_GERG(x, tolerance);
      PseudoCriticalPointGERG(x_GERG,
                              pseudo_critical_point.Tcx,
                              pseudo_critical_point.Dcx);

      if (pseudo_critical_point.Dcx > tolerance)
        pseudo_critical_point.Vcx = 1.0 / pseudo_critical_point.Dcx;

      return pseudo_critical_point;
    }
    // *********************************************************
    template<typename x_type,
             typename value_type>
    auto thermodynamic_properties(const x_type& x,
                                  const gerg_data::Thermodynamic_properties_parameters<value_type>& input_properties,
                                  const value_type& tolerance)
    {
      using size_type = typename x_type::size_type;

      assert(static_cast<size_type>(x.size()) == GERG_num_componets);
      const auto x_GERG = create_x_GERG(x, tolerance);

      gerg_data::Thermodynamic_properties<value_type> thermodynamic_properties;
      thermodynamic_properties.D = value_type();
      thermodynamic_properties.P = value_type();
      thermodynamic_properties.Z = value_type();
      thermodynamic_properties.gamma = value_type();

      int ierr = 0;
      std::string herr;
      DensityGERG(static_cast<int>(input_properties.Type),
                  input_properties.T,
                  input_properties.P,
                  x_GERG,
                  thermodynamic_properties.D,
                  ierr,
                  herr);

      if (ierr != 0)
        throw std::runtime_error(herr);

      PropertiesGERG(input_properties.T,
                     thermodynamic_properties.D,
                     x,
                     thermodynamic_properties.P,
                     thermodynamic_properties.Z,
                     thermodynamic_properties.dPdD,
                     thermodynamic_properties.dPdD2,
                     thermodynamic_properties.d2PdTD,
                     thermodynamic_properties.dPdT,
                     thermodynamic_properties.U,
                     thermodynamic_properties.H,
                     thermodynamic_properties.S,
                     thermodynamic_properties.Cv,
                     thermodynamic_properties.Cp,
                     thermodynamic_properties.W,
                     thermodynamic_properties.G,
                     thermodynamic_properties.JT,
                     thermodynamic_properties.Kappa,
                     thermodynamic_properties.A);

      thermodynamic_properties.gamma = thermodynamic_properties.Cp / thermodynamic_properties.Cv;

      return thermodynamic_properties;
    }
    // *********************************************************
  }
}

#endif
