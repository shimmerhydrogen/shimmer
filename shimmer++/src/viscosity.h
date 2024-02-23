#pragma once
#include "../src/infrastructure_graph.h"


namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


enum viscosity_type
{
    Kukurugya,
    Constant,
};


/// Dynamic viscosity [Pa*s]

    /** Description

    Function for dynamic viscosity calculation 
    of gas with composition X and temperature T.

                            (T1 + 1,47.Tc)       T^1,5
    mu = 10^(-6) . mu1 . X . --------------- . ------------- 
                                T1^1,5        T + 1,47.Tc

    T - Temperature of mixture                                   [K]

    Used coeficients and functions

    TK - constant for conversion Celsius degrees to Kelvins      [K]
    mu1(1).. mu1(13) - dynamic viscosity of gas component 
                        for temperature T1(1)...T(13)            [Pa.s]
    Tc(1) .. Tc(13)  - boiling temperature of gas component      [K]

 
    mu - Dynamic viscosity of gas mixture                        [Pa.s]
          defined for temperature T in [273.15, 2500].
    
    Copyright (C) 2011
    Authors     : Kukurugya Jan, Terpak Jan
    Organization: Technical University of Kosice
    e-mail      : jan.kukurugya@tuke.sk, jan.terpak@tuke.sk
    Revision    : 28.12.2011
    
                  1          2         3        4         5         6        
                 N2         O2       CO2      H2O        SO2       CO        
        mu1 = [  17.580    20.194    14.689    13.022    12.626    17.426    
        T1  = [  20.000    20.000    20.000   120.0      20.0      20.0      
        Tc  = [-195.795  -182.962   -78.400    99.9743  -10.0    -191.51     

				  7          8          9          10         11         12     
				  H2         CH4       C2H4       C2H6       C3H8      n-C4H10  
        mu1 =     8.8129    11.024    10.152      9.2071     8.0115      7.359 
        T1  =    20.000     20.000    20.000     20.000     20.000      20.000   
        Tc  =  -252.754   -161.483  -103.700    -88.598    -42.090      -0.490  
		  
					13         14         15         16        17
				   H2S      i-C4H10    n-C5H12     i-C5H12    C6
        mu1 =     12.674     7.3744     7.2967     7.3340    7.6545 ];
        T1  =     33.000    20.000     40.000     28.000    70.000  ];
        Tc  =    -60.000   -11.749     36.060     27.820    68.710  ];
*/    

double viscosity_kukurugya(const double& temperature,
                 const vector_t& x )
{

    vector_t T1_ratio(17), TcK(17), mu1(17);


    TcK <<   77.355,  90.188, 194.750, 373.1243,  263.15,  81.64,  20.396,
            111.667, 169.450, 184.552, 231.0600,  272.66, 213.15, 261.401,
            309.210, 300.970, 341.860;

    mu1 <<  1.7580e-05,   2.0194e-05,   1.4689e-05,   1.3022e-05,
            1.2626e-05,   1.7426e-05,   8.8129e-06,   1.1024e-05,
            1.0152e-05,   9.2071e-06,   8.0115e-06,   7.3590e-06,
            1.2674e-05,   7.3744e-06,   7.2967e-06,   7.3340e-06,   7.6545e-06;

    /*   Ratio for T1 :  
                        (T1 + 1,47.Tc)   
                        ---------------  
                            T1^1,5       
    */ 
    T1_ratio <<
    8.106103073298680e-02,   8.481949721214364e-02,   1.154431060326530e-01,
    1.207949375457145e-01,   1.354757645828227e-01,   8.231600064362242e-02,
    6.437915648797050e-02,   9.111016202207192e-02,   1.080333653700903e-01,
    1.124563655078777e-01,   1.260774018214931e-01,   1.382610070215963e-01,
    1.156448359952611e-01,   1.349635259891934e-01,   1.385340478072904e-01,
    1.422823416765382e-01,   1.330399925405980e-01;

    //  Check Temperature limits (0°C <  T < 2226.85°C)
    double T = (temperature < 273.15)? 273.15 :
                             (temperature > 2500)? 2500 : temperature;

    /*    
            (T + 1,47.Tc)   
           ---------------  
               T^1,5       
    */ 
    vector_t T_ratio = T * vector_t::Ones(17) + 1.47 * TcK;
    T_ratio /= std::pow(T, 1.5);

    // Warning:  X_factor is equal for each pipe, since temperature is not varying,
    //           So, it could be computed just once.  

    vector_t X_factor = mu1.array() * T1_ratio.array() / T_ratio.array();

    double mu = X_factor(0) * x(GAS_TYPE::N2)  + X_factor(1) * x(GAS_TYPE::O2)  + 
                X_factor(2) * x(GAS_TYPE::CO2) + X_factor(3) * x(GAS_TYPE::H2O) + 
                X_factor(4) * 0.0              + X_factor(5) * x(GAS_TYPE::CO)  + 
                X_factor(6) * x(GAS_TYPE::H2)  + X_factor(7) * x(GAS_TYPE::CH4) + 
                X_factor(8) * 0.0              + X_factor(9) * x(GAS_TYPE::C2H6)+ 
                X_factor(10)* x(GAS_TYPE::C3H8)+ X_factor(11)* x(GAS_TYPE::n_C4H10)+ 
                X_factor(12)* x(GAS_TYPE::H2S) + X_factor(13)* x(GAS_TYPE::i_C4H10)+ 
                X_factor(14)* x(GAS_TYPE::n_C5H12)+ 
                X_factor(15)* x(GAS_TYPE::i_C5H12)+ 
                X_factor(16)*(x(GAS_TYPE::C6H14)  + x(GAS_TYPE::C7H16)  
                                + x(GAS_TYPE::C8H18) + x(GAS_TYPE::C9H20)
                                 + x(GAS_TYPE::C10H22));
    return mu;
}


template<int viscosity_type>
vector_t
viscosity(double temperature, const infrastructure_graph& g)
{
    vector_t mu(num_edges(g));

    auto edge_range = edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    if constexpr (viscosity_type == Kukurugya)
    {
        size_t i = 0;
        for(auto itor = begin; itor != end; itor++,i++ )
        {
            auto pipe = g[*itor]; 
            auto node_in = g[source(*itor, g)];                        
            mu(i) = viscosity_kukurugya(temperature, node_in.gas_mixture);
        }
    } 
    if constexpr (viscosity_type == Constant)
    {
        mu.setConstant(1e-5); 
    }

    return mu;
}

} //end namespace shimmer