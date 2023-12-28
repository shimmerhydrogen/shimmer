#include "MATLAB_GERG_interface.hpp"

namespace GERG
{
  // *********************************************************
  Matlab_interface* Matlab_interface::singleton = nullptr;
  // *********************************************************
  Matlab_interface::Matlab_interface()
  {

  }
  Matlab_interface::~Matlab_interface()
  {
    if (singleton != nullptr)
    {
      delete singleton;
      singleton = nullptr;
    }
  }
  // *********************************************************
  const Matlab_interface& GERG::Matlab_interface::GetInstance()
  {
    if (singleton == nullptr)
      singleton = new Matlab_interface();

    return *singleton;
  }
  // *********************************************************

}
