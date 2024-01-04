#ifndef __GERG_functions_H
#define __GERG_functions_H

namespace GERG
{
  // iFlag=0;
  // [Tr_b,Dr_b] = ReducingParametersGERG(Aplus'*reshape(CC_gas(:,:,1),dimn,21));
  // [Tcx_b,Dcx_b,Vcx_b]=PseudoCriticalPointGERG(Aplus'*reshape(CC_gas(:,:,1),dimn,21),dimb);
  // [Pcheck, Zm, Den]=PropertiesGERG(iFlag, pm(:,1)/1e3, Tb, Aplus'*reshape(CC_gas(:,:,1),dimn,21),dimb,Tr_b,Dr_b,Tcx_b,Dcx_b,Vcx_b);

  template <class matrix_type>
  struct Reducing_parameters final
  {
      matrix_type Tr; ///< reducing temperature (K)
      matrix_type Dr; ///< reducing density (mol/l)
  };

  template <class matrix_type>
  struct Pseudo_critical_point final
  {
      matrix_type Tcx; ///< pseudo critical point temperature (K)
      matrix_type Dcx; ///< pseudo critical point density (mol/l)
      matrix_type Vcx; ///< pseudo critical point volume (l)
  };

  template <class matrix_type>
  struct Thermodynamic_properties final
  {
      matrix_type P1; ///< a pressure (?)
      matrix_type Z; ///< compressibility factor
      matrix_type D; ///< density (mol/l)
      matrix_type gamma; ///< a constant (?)
  };

  /// \brief Calculate reducing variables.
  /// \param x, composition (mole fraction)
  /// \return the reducing parameters struct
  /// \note Only need to call this if the composition has changed.
  template <class matrix_type>
  Reducing_parameters<matrix_type> reducing_parameters(const matrix_type& x);
  /// \brief Calculate a pseudo critical point
  /// as the mole fraction average
  /// of the critical temperatures and critical volumes
  /// \param x, composition (mole fraction)
  /// \param dimn, the composition size
  /// \return the pseudo critical point struct
  template <class matrix_type>
  Pseudo_critical_point<matrix_type> pseudo_critical_point(const matrix_type& x,
                                                           const unsigned int dimn);
  /// \brief Calculate thermodynamic properties as a function of temperature and density.
  /// \param P, pressure (kPa), size dimn x 1
  /// \param T, temperature (K), size dimn x 1
  /// \param x, composition (mole fraction)
  /// \param dimn, the composition size
  /// \param Tr, Dr, the reducing temperature (K) and the reducing density (mol/l), size dimn x 1
  /// \param Tcx, Dcx, Vcx, the pseudo critical point temperature, density and volume, size dimn x 1
  /// \param parameters, other input parameters
  /// \return the thermodynamic properties struct
  template <class matrix_type, class vector_type,
            class parameters_type>
  Thermodynamic_properties<vector_type> 
    thermodynamic_properties(const vector_type& P,
                            const double& T,
                            const matrix_type& x,
                            const unsigned int dimn,
                            const Reducing_parameters<vector_type>& reducing_parameters,
                            const Pseudo_critical_point<vector_type>& pseudo_critical_point, 
                            const parameters_type& parameters);
}

#endif // __GERG_functions_H
