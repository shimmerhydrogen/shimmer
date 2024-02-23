#pragma once

#include <unordered_map>
#include <Eigen/Sparse>
#include "../src/infrastructure_graph.h"

namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


double
viscosity(const double& Temperature,
          const vector_t & X);


double
friction_factor_average(const edge_properties& pipe ,const double & Temperature,
                        const double & flux, const double & mu);


double 
inertia_resistance( const edge_properties& pipe, const double& dt,
                    const double& mean_pressure) ;

double 
friction_resistance(const  edge_properties& pipe,
                    double temperature,
                    double mu,
                    double c2, 
                    double flux);


} //end namespace shimmer
