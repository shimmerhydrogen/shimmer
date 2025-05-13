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

#pragma once


#include <iostream>
#include <Eigen/Sparse>
#include <iomanip>

#include "../src/variable.h"
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
    vector_t area_pipes_;

    incidence inc_;
    const infrastructure_graph& graph_;

public:
    time_solver(const  infrastructure_graph& g,
                double Tm):
                graph_(g), temperature_(Tm)
    {
        inc_ = incidence(g);
        area_pipes_ = area(g);
    }


    void
    pipe_stations_activation(size_t step, const variable& v)
    {
        auto edge_range = boost::edges(graph_);
        auto begin = edge_range.first;
        auto end = edge_range.second;

        for (auto itor = begin; itor != end; itor++)
        {
            auto pipe = graph_[*itor];

            if (pipe.type == pipe_type::PIPE)
                continue;

            auto& st = pipe.pipe_station;

            auto s = boost::source(*itor, graph_);
            auto t = boost::target(*itor, graph_);

            auto source_node = graph_[s].i_snum;
            auto target_node = graph_[t].i_snum;

            st->activate(step, source_node, target_node, v);
        }
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

        pipe_stations_activation(iter, var_time);

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

        std::ofstream ofs("warnings.txt");

        for(size_t it = 1; it <= num_steps; it++, t+=dt)
        {
            ofs.open("warnings.txt", std::ios::app);
            if(!ofs.is_open())
                throw std::runtime_error("ERROR:  warnings file not opened.");
            ofs<<"========================================================"<< std::endl;
            ofs << "Solving at time ...."<< it <<std::endl;
            ofs.close();

            std::cout<<"========================================================"<< std::endl;
            std::cout<<"========================================================"<< std::endl;
            std::cout<<"========================================================"<< std::endl;
            std::cout << "Solving at time ...."<< it <<std::endl;
            std::cout<<"========================================================"<< std::endl;
            std::cout<<"========================================================"<< std::endl;
            std::cout<<"========================================================"<< std::endl;

            pipe_stations_activation(it, var_);

            size_t ic;
            for(ic = 0; ic <= MAX_CONSTRAINT_ITER; ic++)
            {
                ofs.open("warnings.txt", std::ios::app);
                if(!ofs.is_open())
                    throw std::runtime_error("ERROR:  warnings file not opened.");
                ofs << " * Iteration it ..."<<ic<< std::endl;
                ofs.close();

                std::cout<<"****************************************************************"<< std::endl;
                std::cout<<"****************************************************************"<< std::endl;
                std::cout << " Iteration CONSTRAINTS it ..............."<<ic<< " ... at time "<< it << std::endl;
                std::cout<<"****************************************************************"<< std::endl;
                std::cout<<"****************************************************************"<< std::endl;

                EQ_OF_STATE eos;
                eos.compute_molar_mass(y_nodes, y_pipes);

                // To be finish when it is clear how x changes and modifies mu.
                auto mu = viscosity<viscosity_type>(temperature_, graph_);

                linearized_fluid_solver lfs(it, unsteady, tol, dt, temperature_, mu, inc_, graph_);
                lfs.run(area_pipes_, var_guess_, var_, &eos, ic);

                bool pass_constr = lfs.check_constraints(it); 
                bool pass_control =  lfs.check_controls(it);

                if(pass_constr && pass_control)
                {
                    std::cout<< "++++++++++++++++++**** MODIFIED VARIABLE ****++++++++++++++++++++++ " << std::endl;
                    var_ =  lfs.get_variable();

                    //std::cout<< "VARIABLE : \n";
                    //std::cout<<  var_.make_vector() << std::endl;
                    
                    break;
                }
                //std::cout<< "VARIABLE : \n";
                //std::cout<<  lfs.get_variable().make_vector() << std::endl;
            }

            if(ic == MAX_CONSTRAINT_ITER)
                std::cout << "ERROR: FAILURE to apply HARD constraints. Max number of iterations has been reached.";

            var_in_time.row(it) =  var_.make_vector();

            var_guess_ = var_;
        }

        std::ofstream vfs("../build/unit_tests/var_time_DISMA_flux.dat");
        if(!vfs.is_open())
            std::cout<<"WARNING: var_in_time file has not been opened." << std::endl;

        for(size_t i = 0; i < var_in_time.rows(); i++)
        {
            for(size_t j = 0; j < var_in_time.cols(); j++)
                vfs << std::setprecision(16) << var_in_time(i,j) << " ";
            vfs << std::endl;
        }
        vfs.close();

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