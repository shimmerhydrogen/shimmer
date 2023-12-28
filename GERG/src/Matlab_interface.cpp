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
  const Matlab_interface& GERG::Matlab_interface::get_instance()
  {
    if (singleton == nullptr)
    {
      singleton = new Matlab_interface();

      std::string matlab_shimmer_dir(__FILE__);
      matlab_shimmer_dir = matlab_shimmer_dir.substr(0, matlab_shimmer_dir.find_last_of("\\/"));
      matlab_shimmer_dir = matlab_shimmer_dir + "/../../matlab/Mod_SHIMM_v01";

      matlab::data::ArrayFactory factory;

      std::vector<matlab::data::Array> args({
                                              factory.createCharArray(matlab_shimmer_dir)
                                            });

      const auto result = singleton->matlab().feval(u"addpath",
                                                    args);
    }

    return *singleton;
  }
  // *********************************************************
  std::vector<matlab::data::Array> Matlab_interface::reducing_parameters(const matlab::data::TypedArray<double>& x) const
  {
    matlab::engine::MATLABEngine& matlab = Matlab_interface::get_instance().matlab();

    // Call MATLAB function
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
