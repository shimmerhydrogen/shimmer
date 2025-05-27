#ifndef __SHIMMER_GERG_FUNCTIONS_HPP__
#define __SHIMMER_GERG_FUNCTIONS_HPP__

#include "GERG/GERG2008.h"
#include "shimmer_gerg_data.hpp"
#include <stdexcept>
#include <Eigen/Dense>
#include <vector>

namespace shimmer_gerg
{
namespace gerg_functions
{
    static const unsigned int GERG_num_componets = 21;
    // *********************************************************
    template <typename vector_type>
    std::vector<double> 
    create_x_GERG(const vector_type& x,
                  const double& tolerance)
    {
        assert(x.size() == GERG_num_componets);

        std::vector<double> x_GERG(GERG_num_componets + 1);

        x_GERG[0] = 0.0;
        for (unsigned int i = 0; i < GERG_num_componets; ++i)
            x_GERG[i + 1] = std::abs(x[i]) > tolerance ? x[i] : 0.0;

        return x_GERG;
    }
    // *********************************************************
    inline void setup_GERG() noexcept
    {
        SetupGERG();
    }
    // *********************************************************
    void 
    molar_mass(const Eigen::MatrixXd &mole_frac, 
               const double& tolerance,
               Eigen::VectorXd &Mm);
    // *********************************************************
    gerg_data::Reducing_parameters<double> 
    reducing_parameters(const std::vector<double>& x,
                        const double& tolerance);
    // *********************************************************
    gerg_data::Pseudo_critical_point<double> 
    pseudo_critical_point(const std::vector<double>& x,
                          const double& tolerance);
    // *********************************************************
    gerg_data::Thermodynamic_properties<double> 
    thermodynamic_properties(const std::vector<double>& x,
                             const gerg_data::Thermodynamic_properties_parameters<double>& input_properties,
                             const double& tolerance);
    // *********************************************************
    gerg_data::Thermodynamic_properties<Eigen::VectorXd> 
    thermodynamic_properties(const Eigen::VectorXd& temperature,
                             const Eigen::VectorXd& pressure,
                             const Eigen::MatrixXd& x,
                             const gerg_data::Thermodynamic_properties_parameters<double>::Types & type,
                             const double& tolerance);
}
}

#endif
