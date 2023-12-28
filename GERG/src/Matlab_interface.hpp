#ifndef __Matlab_interface_H
#define __Matlab_interface_H

#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

namespace GERG
{
  class Matlab_interface final
  {
    private:
      Matlab_interface();
      ~Matlab_interface();

      static Matlab_interface* singleton;

      std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;

    public:
      Matlab_interface(Matlab_interface &other) = delete;
      void operator=(const Matlab_interface&) = delete;

      static const Matlab_interface& get_instance();

      inline matlab::engine::MATLABEngine& matlab() const
      { return *matlabPtr; }

      template <class matrix_type>
      inline matlab::data::TypedArray<double> matrix_to_matlab(matlab::data::ArrayFactory& factory,
                                                               const matrix_type& matrix) const
      {
        return factory.createArray(
              {
                static_cast<unsigned int>(matrix.rows()),
                static_cast<unsigned int>(matrix.cols())
              },
              matrix.data(),
              matrix.data() + matrix.size());
      }

      template <class matrix_type>
      inline matrix_type matlab_to_matrix(const matlab::data::TypedArray<double>& matrix) const
      {
        matrix_type result;
        const std::vector<unsigned long int> size = matrix.getDimensions();
        result.resize(size.at(0),
                      size.at(1));
        std::copy(matrix.begin(), matrix.end(), result.data());
        return result;
      }

      std::vector<matlab::data::Array> reducing_parameters(const matlab::data::TypedArray<double>& x) const;
  };
}

#endif // __Matlab_interface_H
