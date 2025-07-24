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
#include <fstream>

#include "solver/variable.h"
#include "solver/gas_law.h"
#include "solver/fluid_solver.h"
#include "solver/incidence_matrix.h"
#include "solver/viscosity.h"
#include "solver/boundary.h"
#include "infra/infrastructure.h"
#include <boost/graph/adjacency_list.hpp>

namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


template<typename EQ_OF_STATE, int viscosity_type>
class qt_solver
{

    double temperature_;

    matrix_t y_;

    vector_t rho_msh_;
    variable var_msh_;
    variable var_msh_guess_;
    vector_t area_pipes_;

    matrix_t rho_msh_in_time_;
    matrix_t var_msh_in_time_;

    incidence msh_inc_;
    const infrastructure& infra_;

    EQ_OF_STATE eos;

public:
    qt_solver(const infrastructure& infrain,
              double Tm):
              infra_(infrain), temperature_(Tm)
    {
        msh_inc_ = incidence(infra_.graph);
        area_pipes_ = area(infra_.graph);
    }


    void
    pipe_stations_activation(size_t step, const variable& v)
    {
        auto [ebegin, eend] = boost::edges(infra_.graph);
        for (auto itor = ebegin; itor != eend; itor++)
        {
            auto pipe = infra_.graph[*itor];
            auto s = boost::source(*itor, infra_.graph);
            auto t = boost::target(*itor, infra_.graph);

            if (pipe.type == pipe_type::PIPE)
                continue;

            auto& st = pipe.pipe_station;
            if (!st) {
                std::cerr << "Invalid pointer" << std::endl;
                exit(-1);
            }

            auto source_node = infra_.graph[s].i_snum;
            auto target_node = infra_.graph[t].i_snum;

            st->activate(step, source_node, target_node, v);
        }
    }


    void
    set_initialization( const variable& var)
    {
        var_msh_guess_ = var;
        //Alternatively, rho should be given when setting initialization
        rho_msh_ = vector_t::Ones(boost::num_edges(infra_.graph));
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
        auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

        size_t iter = 0;
        auto var_time = variable(vector_t::Zero(num_vertices(infra_.graph)),
                                 vector_t::Zero(num_edges(infra_.graph)),
                                 vector_t::Zero(num_vertices(infra_.graph)));

        pipe_stations_activation(iter, var_time);

        linearized_fluid_solver lfs(iter, unsteady, tolerance, dt, temperature_, mu, msh_inc_, infra_.graph);
        lfs.run(area_pipes_, var_guess, var_time, &eos);
        var_msh_guess_ = lfs.get_variable();
        rho_msh_ = eos.density(&lfs);
    }

    
    bool
    qt_net_nodes(double dt, size_t at_step, const variable& var_msh,
                 const  matrix_t& y_msh_nodes, vector_t& lhs_nodes, matrix_t& rhs_nodes)
    {
        // Mass conservation for network nodes
        auto inc_mat = msh_inc_.matrix();


        // 1. Pipes Injection/Ejection
        
        // Here, local pipe quantities are needed for network(original) nodes
        auto inc_smat = msh_inc_.matrix();
        for(size_t iN = 0; iN < infra_.num_original_stations; iN++)
        {
            Eigen::SparseVector<double> node_flux = msh_inc_.matrix().row(iN).cwiseProduct(var_msh.flux);       
            
            // Loop by face
            for(Eigen::SparseVector<double>::InnerIterator it(node_flux); it; ++it)
            {
                auto pipe_num = it.col();

                double flux_eject  = std::max(0.0,-it.value());
                double flux_inject = std::max(0.0, it.value());

                // 1.1 Compute ejection : LHS
                lhs_nodes(iN) += flux_eject; 

                // 1.2 Compute injection: RHS

                // 1.2.1 Y@face:  y has to be approx @face (equivalently on the center
                //                of the pipe). Upwind scheme chosen (explicit).  
                auto ed = infra_.p_i2ed[pipe_num]; 

                auto s = boost::source(ed, infra_.graph);
                auto t = boost::target(ed, infra_.graph);

                auto source_num = infra_.graph[s].i_snum;
                auto target_num = infra_.graph[t].i_snum;

                int upw_num; 
                if(source_num == iN )
                    upw_num = target_num;
                else if (target_num == iN)
                    upw_num = source_num;
                else 
                {
                    std::cout << "ERROR: in QT at  network Node" << iN 
                              << ". Loop over pipes arriving to iN."
                              << "Pipe ("<<pipe_num <<"): from "<< source_num
                              << " to "<<target_num <<"  is incosistent with iN. \n.";
                    return SHIMMER_GENERIC_FAILURE;
                }

                //1.2.2 Sum inj
                for (size_t iC = 0; iC < NUM_GASES; iC++)
                {
                    double y_face = y_msh_nodes(upw_num, iC);
                    rhs_nodes(iN,iC) += flux_inject * y_face; 
                }
            }
        }


        /*
        1.2 Problems
        A. I do not know if gas_mixture are molar or mass fractions
        */      

        // 1.3 External Injection/Ejection          
        auto v_range = boost::vertices(infra_.graph);
        for(auto itor = v_range.first; itor != v_range.second; itor++)
        {
            const auto & node_prop = infra_.graph[*itor]; 
            if (not node_prop.node_station) {
                std::string errstr = "Graph node " + std::to_string(*itor) +
                    " points to an invalid station"; 
                throw std::invalid_argument(errstr);
            }

            if (node_prop.type == station_type::FICTITIOUS_JUNCTION)
                continue;

            auto bnd =  node_prop.node_station->boundary();
            auto L = bnd.value();    

            // Injection
            if (L > 0.0)
            {
                for (size_t iC = 0; iC < NUM_GASES; iC++)
                {
                    rhs_nodes(node_prop.i_snum, iC) += L * y_msh_nodes(node_prop.i_snum, iC);
                }
            }    
            // Ejection
            else
                lhs_nodes(node_prop.i_snum) += L;
        
        }
#if 0     
/* if NODE_ACCUMULATES
        // 1.3 Time term
        // V/c2 (dpdt) => I would rather do  V*(d(p/c2)/dt)  

        vector_t press_before = var_in_time_[it-1].pressure.head(num_nodes);
        vector_t press_now  = var_in_time_[it].pressure.head(num_nodes); 
        vector_t c2_now  = c2_nodes.head(num_nodes);  

        vector_t dp = press_now - press_before;
        vector_t phi_vec = phi_vector(dt, c2_now, graph_).array() * dp.array();

        #if 0
        Problems:
        (A) c2:  out of reach here since it depends on the fluid solver, which is out of scope
            => Options: a. save a lfs/eos when iteration are finished, or save c2
        (B) V/c2 (dpdt) => I would rather do  V*(d(p/c2)/dt) thus I would need also c2 in t^{n} 
        (C) global info: num_nodes/num_pipe of the original network are needed.
        #endif
*/
#endif

    }


    void
    qt_refine_nodes(size_t it, double dt)
    {
        for(const auto& pd : infra_.pipe_discretizations)
        {
            auto dtdx =  (dt/pd.dx); 

            // Density pipes
            vector_t rho  = rho_msh_in_time_.row(it-1);  

    #if 0 
            /// vel [m/s] velocity of the gas within pipes => Flux is needed! area_pipes too
            vector_t vel = flux.cwiseQuotient(area_pipes.cwiseProduct(rho));

            // Iterate by global pipe: get velocity local
            vel_plus_half  = inc.matrix_in() * vel_local; 
            vel_minus_half = inc.matrix_out() * vel_local; 
            vel_loc_node   = 0.5* (vel_plus_half + vel_minus_half_);

            // Coefficients
            a_loc_node = 1.0 + dtdx * dtdx * vel_loc_node.array() * ( vel_plus_half.array() - vel_minus_half.array());
            a_plus = dtdx * vel_plus  * (0.5 - dtdx * vel_plus_half);             
            a_minus= dtdx * vel_minus * (0.5 - dtdx * vel_minus_half);             

            // mass fracctions
            Y_loc_nodes_new = a_loc_node * Y_loc_nodes
                                + a_plus * Y_loc_nodes_plus 
                                + a_minus * Y_loc_nodes_minus;  
    #endif
        }
    }


    void
    transmission(size_t it, 
                 double tol,
                 double dt,         // Most of this should be passed in a config struct 
                 size_t MAX_CONSTRAINT_ITER,
                 equation_of_state & eos)
    {
        bool unsteady = true;

        pipe_stations_activation(it, var_msh_);

        size_t ic;
        for(ic = 0; ic <= MAX_CONSTRAINT_ITER; ic++)
        {
            std::cout<<"****************************************************************"<< std::endl;
            std::cout << " Iteration CONSTRAINTS it ..............."<<ic<< " ... at time "<< it << std::endl;
            std::cout<<"****************************************************************"<< std::endl;

            // To be finished when it is clear how x changes and modifies mu.
            auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

            linearized_fluid_solver lfs(it, unsteady, tol, dt, temperature_, mu, msh_inc_, infra_.graph);
            lfs.run(area_pipes_, var_msh_guess_, var_msh_, &eos, ic);

            bool pass_constr = lfs.check_constraints(it); 
            bool pass_control= lfs.check_controls(it);

            if(pass_constr && pass_control)
            {
                std::cout<< "++++++++++++++++++**** MODIFIED VARIABLE ****++++++++++++++++++++++ " << std::endl;
                var_msh_ =  lfs.get_variable();
                rho_msh_ =  eos.density(&lfs);  
                break;
            }
        }

        if(ic == MAX_CONSTRAINT_ITER)
            std::cout << "ERROR: FAILURE to apply HARD constraints. Max number of iterations has been reached.";

        var_msh_in_time_.row(it) =  var_msh_.make_vector();
        rho_msh_in_time_.row(it) =  rho_msh_;

        var_msh_guess_ = var_msh_;
    }


    void
    advance(double dt,
            size_t num_steps,
            double tol,
            const matrix_t& y_msh_nodes,
            const matrix_t& y_msh_pipes)
    {
        //matrix_t y_nodes = make_mass_fraction(num_nodes);
        //matrix_t y_pipes = inc_.matrix_in().transpose() * y_nodes;

        size_t MAX_CONSTRAINT_ITER = 10;

        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 

        var_msh_in_time_ = matrix_t::Zero(num_steps, num_pipes + 2 * num_nodes);
        rho_msh_in_time_ = matrix_t::Zero(num_steps, num_pipes);

        var_msh_ = var_msh_guess_;
        var_msh_in_time_.row(0) =  var_msh_.make_vector();
        rho_msh_in_time_.row(0) =  rho_msh_.transpose();
        
        // save matrices?  y_in_time[it] = y_; 
        
        EQ_OF_STATE eos;

        double t = 0;
        for(size_t it = 1; it < num_steps; it++, t+=dt)
        {
            std::cout<<"========================================================"<< std::endl;
            std::cout << "Solving at time ...."<< it <<std::endl;
            std::cout<<"========================================================"<< std::endl;

            // 1. Update molar masses inside eos
            eos.compute_molar_mass(y_msh_nodes, y_msh_pipes);

            // 2. Fluid solver
            transmission(it, tol, dt, MAX_CONSTRAINT_ITER, eos);

            // 3. [y_nodes, y_pipes] = quality_tracking();
            //matrix_t lhs_nodes = vector_t::Zero(infra_.num_original_stations, NUM_GASES); 
            vector_t lhs_nodes = vector_t::Zero(infra_.num_original_stations); 
            matrix_t rhs_nodes = vector_t::Zero(infra_.num_original_pipes, NUM_GASES); 

            // 3. 1. Continuity at network nodes 
            qt_net_nodes(it, dt, var_msh_, y_msh_nodes);

            // 3. 2. Continuity at discretized nodes 
            qt_refine_nodes( it, dt);

            // 3. 4 Solve Y^n+1            
            matrix_t lhs_inv =  lhs_nodes.cwiseInverse().asDiagonal();

            y_ = lhs_inv.array() * rhs_nodes.array(); 

            // save matrices? y_in_time[it] = y_; 

        }

        return;
    }

    vector_t
    velocity(const pipe_discretization& pd, const variable& var, const vector_t& rho, const vector_t& area) const
    {
        vector_t vel = vector_t::Zero(pd.pipelist.size());
        
        for (const auto & pipe_num : pd.pipelist)
        {
            vel(pipe_num) = var.flux(pipe_num) / (rho(pipe_num) * area(pipe_num));
        }

        return vel;      
    }


    matrix_t 
    velocity_full() const
    {
        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 
        /// vel [m/s] velocity of the gas within pipes.
        matrix_t flux_in_time = var_msh_in_time_.middleCols(num_nodes, num_pipes);
        matrix_t arho_in_time = rho_msh_in_time_;

        for(auto iP = 0; iP < num_pipes; iP++)
            arho_in_time.col(iP) *= area_pipes_(iP);  
        return flux_in_time.array() / arho_in_time.array();       
    }

    vector_t solution() const {return var_msh_.make_vector();}
    vector_t guess() const {return var_msh_guess_.make_vector();}
    matrix_t solution_full() const{ return var_msh_in_time_; }

};

}
