#ifndef __SHIMMER_GERG_DATA_H
#define __SHIMMER_GERG_DATA_H

namespace shimmer_teqp
{
  namespace gerg_data
  {
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
    struct Thermodynamic_properties_parameters final
    {
        enum struct Types
        {
          Gas_phase = 0, ///< strict pressure solver in gas phase without checks (fastest mode, but output state may not be stable single phase)
          Two_phase = 1, ///< to make checks for possible 2-phase state (result may still not be stable single phase, but many unstable states will be identified)
          Liquid_phase = 2 ///< to search for liquid phase (and make the same checks when iFlag=1)
        };

        Types Type; ///< solver type
        matrix_type P; ///< pressure (kPa)
        matrix_type T; ///< temperature (K)
    };

    template <class matrix_type>
    struct Thermodynamic_properties final
    {
        matrix_type P; ///< the resulting pressure (kPa)
        matrix_type Z; ///< compressibility factor
        matrix_type D; ///< density (mol/m^3)
        matrix_type gamma; ///< Cp/Cv
        matrix_type dPdD; ///< First derivative of pressure with respect to density at constant temperature [kPa/(mol/l)]
        matrix_type dPdD2; ///< Second derivative of pressure with respect to density at constant temperature [kPa/(mol/l)^2]
        matrix_type d2PdTD; ///< Second derivative of pressure with respect to temperature and density [kPa/(mol/l)/K]
        matrix_type dPdT; ///< First derivative of pressure with respect to temperature at constant density (kPa/K)
        matrix_type U; ///< Internal energy (J/mol)
        matrix_type H; ///< Enthalpy (J/mol)
        matrix_type S; ///< Entropy [J/(mol-K)]
        matrix_type Cv; ///< Isochoric heat capacity [J/(mol-K)]
        matrix_type Cp; ///< Isobaric heat capacity [J/(mol-K)]
        matrix_type W; ///< Speed of sound (m/s)
        matrix_type G; ///< Gibbs energy (J/mol)
        matrix_type JT; ///< Joule-Thomson coefficient (K/kPa)
        matrix_type Kappa; ///< Isentropic Exponent
        matrix_type A; ///< Helmholtz energy (J/mol)
    };
  }
}

#endif // __SHIMMER_GERG_DATA_H
