
/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#pragma once


#include <iostream>
#include <Eigen/Sparse>

#include "../src/temporal.h"
#include "../src/gas_law.h"
#include "../src/fluid_solver.h"
#include "../src/incidence_matrix.h"


namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; 


template<typename EQ_OF_STATE>
class time_solver
{

    double temperature_;
    variable var_;
    vector_t inlet_nodes_;
    matrix_t flux_ext_;
    vector_t Pset_;
    vector_t area_pipes_;

    incidence inc_;
    const infrastructure_graph& graph_; 

public:
    time_solver(const  infrastructure_graph& g,
                double Tm,
                const  vector_t& Pset,
                const  matrix_t& flux_ext,
                const  vector_t& inlet_nodes): 
                temperature_(Tm), Pset_(Pset), flux_ext_(flux_ext),
                inlet_nodes_(inlet_nodes), graph_(g)
    {
        inc_ = incidence(g);
        area_pipes_ = area(g);
    }


    void
    set_initialization( const variable& var)
    {
        var_ = var; 
    }

    void
    initialization( const variable& var,
                    double dt, 
                    double tolerance, 
                    const matrix_t& y_nodes, 
                    const matrix_t& y_pipes)
    {

        bool unsteady = false;

        var_ = var;
    
        EQ_OF_STATE eos; 
        eos.compute_molar_mass(y_nodes, y_pipes);

        linearized_fluid_solver lfs(unsteady, tolerance, dt, temperature_, inc_, graph_);
        lfs.run(area_pipes_, inlet_nodes_, Pset_(0), flux_ext_.row(0), var_, &eos);
    }
    

    void
    advance(double dt, 
            size_t num_steps, 
            double tolerance, 
            const matrix_t& y_nodes, 
            const matrix_t& y_pipes)
    {
        //matrix_t y_nodes = make_mass_fraction(num_nodes);
        //matrix_t y_pipes = inc_.matrix_in().transpose() * y_nodes;    

        bool unsteady = true;

        double t = 0;
        for(size_t it = 1; it < num_steps; it++, t+=dt)
        {
            std::cout << "Solving at time ...."<<t<< " with iteration..."<< it<< std::endl;
            EQ_OF_STATE eos; 
            eos.compute_molar_mass(y_nodes, y_pipes);

            linearized_fluid_solver lfs(unsteady, tolerance, dt, temperature_,inc_, graph_);
            lfs.run(area_pipes_, inlet_nodes_, Pset_(it), flux_ext_.row(it),  var_, &eos);

            //p_0(:,it)= var_.pressure;
            //G_0(:,it)= var_.flux;
            //L_0(:,it)= var_.L_rate;
            // This maybe should be inside LP(:,it) = dx.*pm(:,k).*AA./cc2b; %kg

            //G_n = G_0(:,ii);
            //p_n = p_0(:,ii);
        }

        std::cout << " * Pressure : ["; //<< std::endl;
        for (int k = 0; k < var_.pressure.size(); ++k)
            std::cout << var_.pressure[k] << ", "  <<  std::endl ;
        std::cout << "]; "<<std::endl;

        std::cout << " * Flux : ["; // << std::endl;
        for (int k = 0; k < var_.flux.size(); ++k)
            std::cout << var_.flux[k] << ", "  <<  std::endl ;
        std::cout << "]; "<<std::endl;

        std::cout << "L_rate = [";// << std::endl;
        for (int k = 0; k < var_.L_rate.size(); ++k)
            std::cout << var_.L_rate[k] << ", " <<  std::endl ;
        std::cout << "]; "<<std::endl;

    }
};

}