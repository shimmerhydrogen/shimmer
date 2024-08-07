#include <iostream>
#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/cpp/deriv_adapter.hpp"
#include "teqp/cpp/derivs.hpp"
#include "teqp/models/GERG/GERG.hpp"

#include <typeinfo>

int main(int argc, char **argv) 
{
  using namespace teqp::cppinterface;
  using namespace teqp::cppinterface::adapter;

  // auto j = R"(
  // {"kind": "GERG2008resid", "model": {"names": ["methane","ethane"]}}
  // )"_json;
  // std::cout<< "model GERG: " << j.dump(2)<< std::endl;

  //auto model = make_model(j);

  const auto& names = teqp::GERG2008::component_names;
  std::vector<std::string> comps = {"carbondioxide", "methane"};
  auto model = teqp::GERG2008::GERG2008ResidualModel(comps);

  const auto molefrac = (Eigen::ArrayXd(2) << 0.2, 0.8).finished();
  const double T = 300;
  const double rho_molar = 3;

  const auto R = model.R(molefrac);
  std::cout<< "R: " << R << std::endl;
  const auto Tc = model.red.get_Tcvec();
  const auto Vc = model.red.get_vcvec();

  return 0;
}
