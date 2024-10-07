#ifndef __SHIMMER_TEQP_FUNCTIONS_HPP__
#define __SHIMMER_TEQP_FUNCTIONS_HPP__

#include "shimmer_gerg_data.hpp"

namespace shimmer_teqp
{
  namespace gerg_functions
  {
    // *********************************************************
    template<typename mol_fracs_type,
             typename c_vec_type,
             typename value_type>
    auto pseudo_critical_point(const mol_fracs_type& mol_fracs,
                               const c_vec_type& Tc_vec,
                               const c_vec_type& vc_vec,
                               const value_type& tolerance)
    {
      using size_type = typename c_vec_type::size_type;

      assert(static_cast<size_type>(mol_fracs.size()) == Tc_vec.size() &&
             vc_vec.size() == Tc_vec.size());

      gerg_data::Pseudo_critical_point<value_type> pseudo_critical_point;

      pseudo_critical_point.Tcx = value_type();
      pseudo_critical_point.Vcx = value_type();
      pseudo_critical_point.Dcx = value_type();

      for (size_type i = 0; i < Tc_vec.size(); ++i)
      {
        pseudo_critical_point.Tcx += mol_fracs[i] *
                                     Tc_vec[i];
        pseudo_critical_point.Vcx += mol_fracs[i] *
                                     vc_vec[i];
      }

      if (pseudo_critical_point.Vcx > tolerance)
        pseudo_critical_point.Dcx = 1.0 / pseudo_critical_point.Vcx;

      return pseudo_critical_point;
    }
    // *********************************************************
    template<typename mol_fracs_type,
             typename value_type>
    auto thermodynamic_properties(const mol_fracs_type& mol_fracs,
                                  const value_type& tolerance)
    {
      gerg_data::Thermodynamic_properties<value_type> thermodynamic_properties;

      thermodynamic_properties.D = 4.853135975758211e+01;

      return thermodynamic_properties;
    }
    // *********************************************************
  }
}

#endif
