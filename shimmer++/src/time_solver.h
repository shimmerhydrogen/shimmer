
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
    vector_t outlet_nodes_;
    matrix_t flux_ext_;
    vector_t Pset_;
    incidence inc_;

public:
    time_solver(const infrastructure_graph& g,
            double Tm,
            const vector_t& Pset,
            const matrix_t& flux_ext,
            const vector_t& inlet_nodes,
            const vector_t& outlet_nodes): 
            temperature_(Tm), Pset_(Pset), flux_ext_(flux_ext),
            inlet_nodes_(inlet_nodes), outlet_nodes_(outlet_nodes)
    {
        inc_ = incidence(g);
    }

    void
    initialization( const vector_t& Pguess,
                    const vector_t& Gguess, 
                    const vector_t& L_rate)
    {
        // Do steady state init here.
        var_ = variable(Pguess, Gguess, L_rate); 
    }
    
    void
    advance(double dt, 
            double num_steps, 
            double tolerance, 
            const infrastructure_graph& graph,
            const matrix_t& y_nodes, 
            const matrix_t& y_pipes)
    {
        //matrix_t y_nodes = make_mass_fraction(num_nodes);
        //matrix_t y_pipes = inc_.matrix_in().transpose() * y_nodes;    

        vector_t area_pipes = area(graph);

        double t = 0;
        for(size_t it = 1; it < num_steps; it++, t+=dt)
        {

            std::cout << "Solving at time ...."<<t<< " with iteration..."<< it<< std::endl;
            EQ_OF_STATE eos; 
            eos.compute_molar_mass(y_nodes, y_pipes);

            linearized_fluid_solver lfs(tolerance, dt, temperature_,inc_, graph);
            lfs.run(area_pipes, inlet_nodes_, Pset_.row(it), flux_ext_.col(it),  var_, &eos);

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