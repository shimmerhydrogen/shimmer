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
      singleton = new Matlab_interface();

    return *singleton;
  }
  // *********************************************************
  void Matlab_interface::reducing_parameters(const matlab::data::TypedArray<double>& x) const
  {
    matlab::engine::MATLABEngine& matlab = Matlab_interface::get_instance().matlab();

    // Call MATLAB function
    const auto results = matlab.feval(u"ReducingParametersGERG", x);
  }
  // *********************************************************
}
