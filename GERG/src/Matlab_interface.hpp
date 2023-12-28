#ifndef __Matlab_interface_H
#define __Matlab_interface_H

#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include <list>

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

      static std::string matlab_shimmer_directory_path()
      {
        std::string matlab_shimmer_dir(__FILE__);
        matlab_shimmer_dir = matlab_shimmer_dir.substr(0, matlab_shimmer_dir.find_last_of("\\/"));
        matlab_shimmer_dir = matlab_shimmer_dir + "/../../matlab/Mod_SHIMM_v01";
        return matlab_shimmer_dir;
      }

      static const Matlab_interface& get_instance(const std::list<std::string>& add_directory_paths = {});

      inline matlab::engine::MATLABEngine& engine() const
      { return *matlabPtr; }

      template <class matrix_type>
      static matlab::data::TypedArray<double> matrix_to_matlab(matlab::data::ArrayFactory& factory,
                                                               const matrix_type& matrix)
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
      static matrix_type matlab_to_matrix(const matlab::data::TypedArray<double>& matrix)
      {
        matrix_type result;
        const std::vector<unsigned long int> size = matrix.getDimensions();
        result.resize(size.at(0),
                      size.at(1));
        std::copy(matrix.begin(), matrix.end(), result.data());
        return result;
      }

      static std::u16string string_to_matlab(const std::string& s)
      {
        std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> conv;
        return conv.from_bytes(s);
      }

      std::string matlab_to_string(const std::u16string& s) const
      {
        std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> conv;
        return conv.to_bytes(s);
      }

      bool is_directory_on_matlab_path(const std::string& directory_path) const;

      std::vector<matlab::data::Array> reducing_parameters(const matlab::data::TypedArray<double>& x) const;
  };
}

#endif // __Matlab_interface_H
