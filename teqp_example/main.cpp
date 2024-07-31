#include <iostream>
#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/cpp/deriv_adapter.hpp"
#include "teqp/cpp/derivs.hpp"

int main(int argc, char **argv) 
{
	using namespace teqp::cppinterface;
	using namespace teqp::cppinterface::adapter;

	auto j = R"(
	{"kind": "GERG2008resid", "model": {"names": ["methane","ethane"]}}
	)"_json;
	std::cout<< "j: " << j.dump(2)<< std::endl;
	auto ptr = make_model(j);

	return 0;
}
