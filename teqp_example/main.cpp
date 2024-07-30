#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/cpp/deriv_adapter.hpp"

int main(int argc, char **argv) 
{
	using namespace teqp::cppinterface;
	using namespace teqp::cppinterface::adapter;

	//# Note that names are case-sensitive; this doesn't work
	//model = teqp.make_model({'kind':"GERG2008resid", 'model':{"names": ['methane','ethane']}})

	auto j = R"(
	{"kind": "GERG2008resid", "model": {"names": ['methane','ethane']}}
	)"_json;
	auto ptr = make_model(j);

	return 0;
}
