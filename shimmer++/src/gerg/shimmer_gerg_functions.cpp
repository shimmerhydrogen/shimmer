#include "shimmer_gerg_functions.hpp"
#include <iostream>
namespace shimmer_gerg
{
namespace gerg_functions
{
// *********************************************************
void molar_mass(const Eigen::MatrixXd &x,                 
                const double& tolerance,
                Eigen::VectorXd &Mm)
{
    Mm = Eigen::VectorXd::Zero(x.rows());
    for(auto iR = 0; iR < x.rows(); iR++)
    {
        const auto x_GERG = create_x_GERG(x.row(iR), tolerance);        
        MolarMassGERG(x_GERG, Mm(iR));
    }
}

// *********************************************************
gerg_data::Reducing_parameters<double> 
reducing_parameters(const std::vector<double>& x,
                    const double& tolerance)
{
    assert(x.size() == GERG_num_componets);

    const auto x_GERG = create_x_GERG(x, tolerance);
    gerg_data::Reducing_parameters<double> reducing_parameters;

    reducing_parameters.Tr = double();
    reducing_parameters.Dr = double();

    ReducingParametersGERG(x_GERG,
                            reducing_parameters.Tr,
                            reducing_parameters.Dr);

    return reducing_parameters;
}
// *********************************************************
gerg_data::Pseudo_critical_point<double> 
pseudo_critical_point( const std::vector<double>& x,
                       const double& tolerance)
{
    assert(x.size() == GERG_num_componets);

    gerg_data::Pseudo_critical_point<double> pseudo_critical_point;

    pseudo_critical_point.Tcx = double();
    pseudo_critical_point.Vcx = double();
    pseudo_critical_point.Dcx = double();

    const auto x_GERG = create_x_GERG(x, tolerance);

    PseudoCriticalPointGERG(x_GERG,
                            pseudo_critical_point.Tcx,
                            pseudo_critical_point.Dcx);

    if (pseudo_critical_point.Dcx > tolerance)
        pseudo_critical_point.Vcx = 1.0 / pseudo_critical_point.Dcx;

    return pseudo_critical_point;
}
// *********************************************************
gerg_data::Thermodynamic_properties<double>
thermodynamic_properties(const std::vector<double>& x,
                         const gerg_data::Thermodynamic_properties_parameters<double>& input_properties,
                         const double& tolerance)
{
    assert(x.size() == GERG_num_componets);
    const auto x_GERG = create_x_GERG(x, tolerance);

    // Shimmer++ uses Internation System Units, so pressure is given in [Pa].
    // As AGA8CODE uses [kPa] a rescaling must be done  
    double PressureScaling = 1.e-03;

    gerg_data::Thermodynamic_properties<double> thermodynamic_properties;
    thermodynamic_properties.D = double();
    thermodynamic_properties.P = double();
    thermodynamic_properties.Z = double();
    thermodynamic_properties.gamma = double();

    int ierr = 0;
    std::string herr;
    DensityGERG(static_cast<int>(input_properties.Type),
                input_properties.T,
                input_properties.P * PressureScaling,
                x_GERG,
                thermodynamic_properties.D,
                ierr,
                herr);

    if (ierr != 0)
        throw std::runtime_error(herr);

    PropertiesGERG(input_properties.T,
                    thermodynamic_properties.D,
                    x_GERG,
                    thermodynamic_properties.P,
                    thermodynamic_properties.Z,
                    thermodynamic_properties.dPdD,
                    thermodynamic_properties.dPdD2,
                    thermodynamic_properties.d2PdTD,
                    thermodynamic_properties.dPdT,
                    thermodynamic_properties.U,
                    thermodynamic_properties.H,
                    thermodynamic_properties.S,
                    thermodynamic_properties.Cv,
                    thermodynamic_properties.Cp,
                    thermodynamic_properties.W,
                    thermodynamic_properties.G,
                    thermodynamic_properties.JT,
                    thermodynamic_properties.Kappa,
                    thermodynamic_properties.A);

    thermodynamic_properties.gamma = thermodynamic_properties.Cp / thermodynamic_properties.Cv;

    return thermodynamic_properties;
}
// *********************************************************
gerg_data::Thermodynamic_properties<Eigen::VectorXd> 
thermodynamic_properties(const Eigen::VectorXd& temperature,
                         const Eigen::VectorXd& pressure,
                         const Eigen::MatrixXd& x,
                         const gerg_data::Thermodynamic_properties_parameters<double>::Types & type,
                         const double& tolerance)
{
    assert(x.cols() == GERG_num_componets);

    // Shimmer++ uses Internation System Units, so pressure is given in [Pa].
    // As AGA8CODE uses [kPa] a rescaling must be done  
    double PressureScaling = 1.0e-3;
    
    gerg_data::Thermodynamic_properties<Eigen::VectorXd> thermodynamic_properties;

    auto rows = x.rows();

    thermodynamic_properties.D = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.P = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.Z = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.gamma = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.dPdD= Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.dPdD2 = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.d2PdTD = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.dPdT = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.U  = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.H  = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.S  = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.Cv = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.Cp = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.W  = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.G  = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.JT = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.Kappa = Eigen::VectorXd::Zero(rows);
    thermodynamic_properties.A  = Eigen::VectorXd::Zero(rows);
                        
    for(auto iR = 0; iR <  x.rows(); iR++)
    {
        const auto x_GERG = create_x_GERG(x.row(iR), tolerance);
        /*
        std::cout << "x_row: "<< x.row(iR) << std::endl;
        std::cout << "temp : "<< temperature[iR] << std::endl;
        std::cout << "press: "<< pressure[iR]*PressureScaling<< " [kPa]" << std::endl;
        std::cout << "type: "<< static_cast<int>(type) << std::endl;
        std::cout << "X_gerg: ";
        for(const auto& xv : x_GERG) 
            std::cout << " "<< xv; 
        std::cout<<std::endl;
        std::cout << "tol : "<< tolerance << std::endl;
        */
        int ierr = 0;
        std::string herr;
        DensityGERG(static_cast<int>(type),
                    temperature[iR],
                    pressure[iR]*PressureScaling,
                    x_GERG,
                    thermodynamic_properties.D[iR],
                    ierr,
                    herr);

        if (ierr != 0)
            throw std::runtime_error(herr);

        PropertiesGERG( temperature[iR],
                        thermodynamic_properties.D[iR],
                        x_GERG,
                        thermodynamic_properties.P[iR],
                        thermodynamic_properties.Z[iR],
                        thermodynamic_properties.dPdD[iR],
                        thermodynamic_properties.dPdD2[iR],
                        thermodynamic_properties.d2PdTD[iR],
                        thermodynamic_properties.dPdT[iR],
                        thermodynamic_properties.U[iR],
                        thermodynamic_properties.H[iR],
                        thermodynamic_properties.S[iR],
                        thermodynamic_properties.Cv[iR],
                        thermodynamic_properties.Cp[iR],
                        thermodynamic_properties.W[iR],
                        thermodynamic_properties.G[iR],
                        thermodynamic_properties.JT[iR],
                        thermodynamic_properties.Kappa[iR],
                        thermodynamic_properties.A[iR]);

        thermodynamic_properties.gamma[iR] = thermodynamic_properties.Cp[iR] / thermodynamic_properties.Cv[iR];

        //std::cout << "Z : "<<  thermodynamic_properties.Z[iR] << std::endl;
    }
    return thermodynamic_properties;
}

}
}