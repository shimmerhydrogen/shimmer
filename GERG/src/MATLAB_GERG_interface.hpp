#ifndef __MATLAB_GERG_interface_H
#define __MATLAB_GERG_interface_H

#include "GERG_functions.hpp"

namespace GERG
{
  class Matlab_interface final
  {
    private:
      Matlab_interface();
      ~Matlab_interface();

      static Matlab_interface* singleton;

    public:
      Matlab_interface(Matlab_interface &other) = delete;
      void operator=(const Matlab_interface&) = delete;

      static const Matlab_interface& GetInstance();

      template <class matrix_type>
      Reducing_parameters<matrix_type> reducing_parameters(const matrix_type& x) const
      {
        Reducing_parameters<matrix_type> result;

        return result;
      }
  };
}

#endif // __MATLAB_GERG_interface_H
