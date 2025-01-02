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

        Types Type;
        matrix_type P;
        matrix_type T;
    };

    template <class matrix_type>
    struct Thermodynamic_properties final
    {
        matrix_type P1; ///< a pressure (?)
        matrix_type Z; ///< compressibility factor
        matrix_type D; ///< density (mol/m^3)
        matrix_type gamma; ///< a constant (?)
    };
  }
}

#endif // __SHIMMER_GERG_DATA_H
