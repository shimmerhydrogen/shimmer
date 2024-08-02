#include <iostream>
#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/cpp/deriv_adapter.hpp"
#include "teqp/cpp/derivs.hpp"

#include <typeinfo>

int main(int argc, char **argv) 
{
  using namespace teqp::cppinterface;
  using namespace teqp::cppinterface::adapter;

  auto j = R"(
	{"kind": "GERG2008resid", "model": {"names": ["methane","ethane"]}}
	)"_json;
  std::cout<< "model GERG: " << j.dump(2)<< std::endl;

  auto model = make_model(j);

  const auto molefrac = (Eigen::ArrayXd(2) << 0.2, 0.8).finished();
  const double T = 300;
  const double rho_molar = 3;

  const auto R = model->R(molefrac);
  std::cout<< "R: " << R << std::endl;
  const auto Ar_01 = model->get_Arxy(0, 1, T, rho_molar, molefrac);
  std::cout<< "Ar_01: " << Ar_01 << std::endl;

  return 0;
}
