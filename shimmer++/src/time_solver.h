
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
#include <iomanip>

#include "../src/temporal.h"
#include "../src/gas_law.h"
#include "../src/fluid_solver.h"
#include "../src/incidence_matrix.h"
#include "../src/viscosity.h"
#include "../src/boundary.h"

namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; 


template<typename EQ_OF_STATE, int viscosity_type>
class time_solver
{

    double temperature_;
    variable var_guess_;
    variable var_;
    matrix_t flux_ext_;
    vector_t area_pipes_;

    incidence inc_;
    const infrastructure_graph& graph_; 

public:
    time_solver(const  infrastructure_graph& g,
                double Tm,
                const  matrix_t& flux_ext): 
                graph_(g), temperature_(Tm), 
                flux_ext_(flux_ext)
    {
        inc_ = incidence(g);
        area_pipes_ = area(g);    
    }


    void
    set_initialization( const variable& var)
    {
        var_guess_ = var; 
    }



    void
    initialization( const variable& var_guess,
                    double dt, 
                    double tolerance, 
                    const matrix_t& y_nodes, 
                    const matrix_t& y_pipes)
    {

        bool unsteady = false;

    
        EQ_OF_STATE eos; 
        eos.compute_molar_mass(y_nodes, y_pipes);

        // To be finish when it is clear how x changes and modifies mu. 
        auto mu = viscosity<viscosity_type>(temperature_, graph_); 

        size_t iter = 0;
        auto var_time = variable(vector_t::Zero(num_vertices(graph_)),
                                 vector_t::Zero(num_edges(graph_)),
                                 vector_t::Zero(num_vertices(graph_)));
        linearized_fluid_solver lfs(iter, unsteady, tolerance, dt, temperature_, mu, inc_, graph_);
        lfs.run(area_pipes_, var_guess, var_time, &eos);
        var_guess_ = lfs.get_variable(); 
    }




    void
    advance(double dt, 
            size_t num_steps, 
            double tol, 
            const matrix_t& y_nodes, 
            const matrix_t& y_pipes)
    {
        //matrix_t y_nodes = make_mass_fraction(num_nodes);
        //matrix_t y_pipes = inc_.matrix_in().transpose() * y_nodes;    

        size_t MAX_CONSTRAINT_ITER = 10;

        bool unsteady = true;

        matrix_t var_in_time(num_steps +1, num_edges(graph_) + 2 * num_vertices(graph_));
 
        var_ = var_guess_;
        var_in_time.row(0) =  var_.make_vector();

        double t = 0;
        for(size_t it = 1; it < num_steps; it++, t+=dt)
        {  
            std::ofstream ofs;
            ofs.open("warnings.txt", std::ios::app);
            ofs<<"========================================================"<< std::endl;
            ofs << "Solving at time ...."<< it <<std::endl;
            ofs.close();

            std::cout<<"========================================================"<< std::endl;
            std::cout << "Solving at time ...."<< it <<std::endl;
            std::cout<<"========================================================"<< std::endl;

            size_t ic; 
            for(ic = 0; ic <= MAX_CONSTRAINT_ITER; ic++)
            {
                std::ofstream ofs;
                ofs.open("warnings.txt", std::ios::app);
                ofs << " * Iteration it ..."<<ic<< std::endl;
                ofs.close();

                std::cout<<"***************************************************"<< std::endl;
                std::cout << " Iteration it ..."<<ic<< " ... at time "<< it << std::endl;
                std::cout<<"***************************************************"<< std::endl;

                EQ_OF_STATE eos; 
                eos.compute_molar_mass(y_nodes, y_pipes);

                // To be finish when it is clear how x changes and modifies mu. 
                auto mu = viscosity<viscosity_type>(temperature_, graph_); 

                linearized_fluid_solver lfs(it, unsteady, tol, dt, temperature_, mu, inc_, graph_);
                lfs.run(area_pipes_, var_guess_, var_, &eos);  
                if(lfs.check_constraints(it))
                {
                    std::cout<< "++++++++++++++++++**** MODIFIED VARIABLE ****++++++++++++++++++++++ " << std::endl;
                    var_ =  lfs.get_variable();
                    break;
                }
            }
            
            if(ic == MAX_CONSTRAINT_ITER)
                std::cout << "ERROR: FAILURE to apply HARD constraints. Max number of iterations has been reached.";

            var_in_time.row(it) =  var_.make_vector();

            //var_guess_ = var_;
        }

        std::ofstream ofs("var_in_time.dat");
        if(!ofs.is_open())
            std::cout<<"WARNING: var_in_time file has not been opened." << std::endl;

        for(size_t i = 0; i < var_in_time.rows(); i++)
        {
            for(size_t j = 0; j < var_in_time.cols(); j++)
                ofs << std::setprecision(16) << var_in_time(i,j) << " ";
            ofs << std::endl; 
        }
        ofs.close();

/*
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
*/
        return;    
    }

    vector_t solution(){return var_.make_vector();}
    vector_t guess(){return var_guess_.make_vector();}
};

}