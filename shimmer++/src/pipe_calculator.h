#pragma once

#include <unordered_map>
#include <Eigen/Sparse>
#include "../src/infrastructure_graph.h"

namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


double
viscosity(const double& Temperature,
          const std::unordered_map<std::string, double> & X);

vector_t
friction_factor_average(const double & Temperature, const vector_t & flux, 
                        const infrastructure_graph & graph);
        
} //end namespace shimmer
