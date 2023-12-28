#include "Matlab_interface.hpp"

namespace GERG
{
  // *********************************************************
  Matlab_interface* Matlab_interface::singleton = nullptr;
  // *********************************************************
  Matlab_interface::Matlab_interface() :
    matlabPtr(matlab::engine::startMATLAB())
  {
  }
  Matlab_interface::~Matlab_interface()
  {
  }
  // *********************************************************
  const Matlab_interface& GERG::Matlab_interface::get_instance(const std::list<std::string>& add_directory_paths)
  {
    if (singleton == nullptr)
    {
      singleton = new Matlab_interface();

      matlab::data::ArrayFactory factory;

      for (const std::string& add_directory_path : add_directory_paths)
      {
        std::vector<matlab::data::Array> args({
                                                factory.createCharArray(add_directory_path)
                                              });

        singleton->engine().feval(u"addpath",
                                  args);
      }
    }

    return *singleton;
  }
  // *********************************************************
  bool Matlab_interface::is_directory_in_matlab_path(const std::string& directory_path) const
  {
    matlab::engine::MATLABEngine& matlab = Matlab_interface::get_instance().engine();

    std::string matlab_command =
        "hasFolder = ~isempty(strfind(path, ['" + directory_path + "', pathsep]));";

    matlab.eval(string_to_matlab(matlab_command));
    matlab::data::TypedArray<bool> hasFolder = matlab.getVariable(u"hasFolder");

    return hasFolder[0];
  }
  // *********************************************************
  std::vector<matlab::data::Array> Matlab_interface::reducing_parameters(const matlab::data::TypedArray<double>& x) const
  {
    matlab::engine::MATLABEngine& matlab = Matlab_interface::get_instance().engine();

    std::vector<matlab::data::Array> args({
                                            x
                                          });

    const unsigned int numReturned = 2;
    const std::vector<matlab::data::Array> results = matlab.feval(u"ReducingParametersGERG",
                                                                  numReturned,
                                                                  args);

    return results;
  }
  // *********************************************************
}
