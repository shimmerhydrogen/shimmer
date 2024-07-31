#include <iostream>
#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/cpp/deriv_adapter.hpp"
#include "teqp/cpp/derivs.hpp"

int main(int argc, char **argv) 
{
	using namespace teqp::cppinterface;
	using namespace teqp::cppinterface::adapter;

	//# Note that names are case-sensitive; this doesn't work
	//model = teqp.make_model({'kind':"GERG2008resid", 'model':{"names": ['methane','ethane']}})

	//auto j = R"(
	//{"kind": "GERG2008resid", "model": {"names": ['methane','ethane']}}
	//)"_json;
	
	/*nlohmann::json j = { 
        {"kind", "GERG2008resid"}, 
        {"model", {"names", ["methane","ethane"]}}
    };
    std::cout << j.dump(2);
	auto ptr = make_model(j);*/
	
	nlohmann::json j = { 
        {"kind", "multifluid"}, 
        {"model", {
            {"components", {"../mycp/dev/fluids/Methane.json","../mycp/dev/fluids/Ethane.json"}},
            {"BIP", "../mycp/dev/mixtures/mixture_binary_pairs.json"},
            {"departure", "../mycp/dev/mixtures/mixture_departure_functions.json"}
        }
    }};
    //std::cout << j.dump(2);
    auto am = teqp::cppinterface::make_model(j);

	return 0;
}
