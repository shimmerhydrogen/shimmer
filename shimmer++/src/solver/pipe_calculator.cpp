/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <iostream>

#include "solver/pipe_calculator.h" 
namespace shimmer{

/* 

vector_t
friction_factor_average(const double & Temperature, const vector_t & flux, 
                        const infrastructure_graph & graph)
{
    
    vector_t f(num_edges(graph));

    size_t i = 0;
    auto edge_range = edges(graph);
    auto begin = edge_range.first;
    auto end = edge_range.second;
    for(auto itor = begin; itor != end; itor++,i++ )
    {
        auto node_in = source(*itor, graph);
        auto mu = viscosity(Temperature, graph[node_in].gas_mixture);
        auto pipe = graph[*itor];   
        auto Re = std::abs(flux(i)) * pipe.diameter / (pipe.area() * mu) ; 

        auto eps_over_d = pipe.friction_factor / pipe.diameter;
        auto a = 1.0 / (1.0 + std::pow( Re  / 2720.0, 9));
        auto b = 1.0 / (1.0 +  std::pow(Re * eps_over_d/160.0, 2.0) );

        auto t0 = 3.7 / eps_over_d;
        auto t1 = std::pow( 64.0 / Re, a) ;
        auto t2 = std::pow( 1.8 * std::log10(Re/6.8),  2.0 * (a -1.0) * b);
        auto t3 = std::pow( 2.0 * std::log10(t0), 2.0 * (a - 1.0) * (1.0 - b));
        
        f(i) = t1 * t2 * t3;
    }

    return f; 
}
*/

/// Friction factor obtain from the numerical solution of the Colebrook-White equation
double
friction_factor_average(const edge_properties& pipe, const double & Temperature,
                        const double & flux, const double & mu)
{
    // Friction average cannot be computed when flux = 0, since Re = 0 
    // Therefore a correction is done in impose edge model (fluid solver)
    // Using units of [Kg/s], we consider a flux equal zero when equal to 1.e-12
    if(std::abs(flux) <= 1.e-12) 
        return 0.;

    auto Re = std::abs(flux) * pipe.diameter / (pipe.area() * mu) ; 

    auto eps_over_d = pipe.friction_factor / pipe.diameter;
    auto a = 1.0 / (1.0 + std::pow( Re  / 2720.0, 9));
    auto b = 1.0 / (1.0 +  std::pow(Re * eps_over_d/160.0, 2.0) );

    auto t0 = 3.7 / eps_over_d;
    auto t1 = std::pow( 64.0 / Re, a) ;
    auto t2 = std::pow( 1.8 * std::log10(Re/6.8),  2.0 * (a -1.0) * b);
    auto t3 = std::pow( 2.0 * std::log10(t0), 2.0 * (a - 1.0) * (1.0 - b));
        
    return t1 * t2 * t3; 
}

double 
inertia_resistance( const edge_properties& pipe, const double& dt,
                    const double& mean_pressure) 
{
    return  2.0 * pipe.length * mean_pressure / (dt * pipe.area()); 
}


double 
friction_resistance(const  edge_properties& pipe,
                    double temperature,
                    double mu,
                    double c2, 
                    double flux) 
{
    auto f = friction_factor_average(pipe, temperature, flux, mu);
    auto a = pipe.area();

    //std::cout << f  << "*" << c2 << "*" << pipe.length << "/ ("<< a<< "*" <<a<< "*" <<pipe.diameter <<")" << std::endl;

    return  f * c2 * pipe.length / (a * a * pipe.diameter );    
}




} //end namespace shimmer
