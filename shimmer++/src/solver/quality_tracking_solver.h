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
using vector_row_t = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


template<typename EQ_OF_STATE, int viscosity_type>
class qt_solver
{
    bool refine_;
    double temperature_;
    int num_time_steps_; 
    int MAX_ITERS_STEADY_;
    int TOL_MASSFRAC_;

    matrix_t massfrac_guess_;

    vector_t rho_all_pipes_;
    vector_t rho_all_nodes_;
    variable var_all_;
    variable var_all_guess_;
    vector_t area_all_pipes_;

    matrix_t rho_all_pipes_evol_;
    matrix_t rho_all_nodes_evol_;
    matrix_t var_all_evol_;
    std::vector<matrix_t> molfrac_all_evol_;

    incidence inc_all_;
    infrastructure& infra_;

public:
    qt_solver(infrastructure& infrain,
              double Tm,
              const bool& refine, 
              int num_steps):
              infra_(infrain), temperature_(Tm), refine_(refine), num_time_steps_(num_steps)
    {
        MAX_ITERS_STEADY_ = 500;
        TOL_MASSFRAC_ = 1.e-4; 
        inc_all_ = incidence(infra_.graph);
        area_all_pipes_ = area(infra_.graph);

        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 

        var_all_evol_ = matrix_t::Zero(num_steps, num_pipes + 2 * num_nodes);
        rho_all_pipes_evol_ = matrix_t::Ones(num_steps, num_pipes);
        rho_all_nodes_evol_ = matrix_t::Ones(num_steps, num_nodes);
        molfrac_all_evol_ = std::vector<matrix_t>(num_steps);
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


    template<typename Derived>
    void 
    clip_and_renormalize(Eigen::MatrixBase<Derived>&& mat)
    {
        if(mat.cols() != NUM_GASES)
            throw std::invalid_argument("ERROR:QT: block to clip and renormalized must have #cols = NUM_GASES.");
/*
        mat.unaryExpr([&](double val){

            if(std::abs(val) < TOL_MASSFRAC_)  
            {
                exit(2); 
                return 0.;
            }
            else if(val < 0.) 
            {
                std::cerr << "ERROR:QT: Negative mass fractions" << std::endl;
                exit(1);  
            }
            return val;
        });
*/
        for (int iRow = 0; iRow < mat.rows(); iRow++)
        {
            for (int iCol = 0; iCol < mat.cols(); iCol++)
            {
                if(std::abs(mat(iRow, iCol)) < TOL_MASSFRAC_)  
                    mat(iRow, iCol) = 0; 
                else if(mat(iRow, iCol) < 0.) 
                {
                    std::cerr<< mat(iRow, iCol) << std::endl;
                    throw std::invalid_argument("ERROR:QT: Negative mass fractions"); 
                }
                else  
                    continue; 
            }
        }


        //Renormalize
        for (int iRow = 0; iRow < mat.rows(); iRow++)
        {
            double sum = (mat.row(iRow)).sum();            

            if(std::abs(sum) < TOL_MASSFRAC_)
                throw std::invalid_argument("ERROR:QT: mass fractions sum is close to zero ");                

            mat.row(iRow) /= sum;
        }

        return;
    }


    void
    initialization( const variable& var_guess,
                    double dt,
                    double tolerance)
    {
        bool unsteady = false;

        EQ_OF_STATE eos;    
        auto [massfrac_nodes_, massfrac_pipes_] =  eos.molfrac_2_massfrac(infra_.graph, inc_all_);
        // Mass fraction by comp, needs total molar mass. So molar mass has to be updated any time molfrac changes!

        // Viscosity changes with molfrac (stored by comp in the graph).
        auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

        size_t iter = 0;
        auto var_time = variable(vector_t::Zero(num_vertices(infra_.graph)),
                                 vector_t::Zero(num_edges(infra_.graph)),
                                 vector_t::Zero(num_vertices(infra_.graph)));

        pipe_stations_activation(iter, var_time);

        linearized_fluid_solver lfs(iter, unsteady, tolerance, dt, temperature_, mu, inc_all_, infra_.graph);
        lfs.run(area_all_pipes_, var_guess, var_time, &eos);
        var_all_guess_ = lfs.get_variable();
        var_all_   = var_all_guess_;
        rho_all_pipes_  = eos.density(&lfs);
        rho_all_nodes_  = eos.density_nodes(&lfs);
        massfrac_guess_ = massfrac_nodes_;

        var_all_evol_.row(0) =  var_all_.make_vector();
        rho_all_pipes_evol_.row(0) =  rho_all_pipes_.transpose();
        rho_all_nodes_evol_.row(0) =  rho_all_nodes_.transpose();
    }

    
    vector_row_t 
    flux_exact(const vector_row_t& TransportedVariable, const double& Velocity)
    {
        return Velocity * TransportedVariable;  
    }

    std::pair<vector_row_t, vector_row_t>
    LaxWendroff(double dtdx, const matrix_t& Q, const matrix_t&  flux_nodes, size_t i)
    {
        vector_row_t QL = 0.5 * (Q.row(i-1) + Q.row(i))  - 0.5 * dtdx * (flux_nodes.row(i) - flux_nodes.row(i-1)); 
        vector_row_t QR = 0.5 * (Q.row(i+1) + Q.row(i))  - 0.5 * dtdx * (flux_nodes.row(i+1) - flux_nodes.row(i));  

        return std::make_pair(QL, QR);
    }


    bool
    admixing(bool unsteady, 
            size_t at_step, 
            double dt,  
            const variable& var_all,
            const matrix_t& massfrac_now,
            matrix_t& massfrac_next)
    {
        int num_nodes = (unsteady)? infra_.num_original_stations : boost::num_vertices(infra_.graph);

        vector_t lhs_nodes = vector_t::Zero(num_nodes);
        matrix_t rhs_nodes = matrix_t::Zero(num_nodes, NUM_GASES );

        // Mass conservation for network nodes
        auto inc_mat = inc_all_.matrix();


        // 1. Pipes Injection/Ejection
        
        using svec_itor_t = Eigen::SparseVector<double>::InnerIterator;

        // Here, local pipe quantities are needed for network(original) nodes
        auto inc_smat = inc_all_.matrix();
        for(size_t iN = 0; iN < num_nodes; iN++)
        {
            Eigen::SparseVector<double> inc_node_i = inc_all_.matrix().row(iN);
            
            // Loop by face
            for(svec_itor_t it(inc_node_i); it; ++it)
            {
                auto  pipe_num = it.row();

                double inc_val = it.value();
                double val = inc_val * var_all.flux(pipe_num); 

                double flux_eject  =  std::max(0.0, val);
                double flux_inject =  std::max(0.0, -val);

                // 1.1 Compute ejection : LHS
                lhs_nodes(iN) += flux_eject; 

                // 1.2 Compute injection: RHS

                // 1.2.1 Y@face:  y has to be approx @face (equivalently on the center
                //                of the pipe). Upwind scheme chosen (explicit).  
                edge_descriptor ed = infra_.p_i2ed[pipe_num]; 

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
                    std::cout << "ERROR: in QT at  network node iN = " << iN 
                              << ". Loop over pipes arriving to this node."
                              << " Pipe ("<< pipe_num <<"): from "<< source_num
                              << " to "<< target_num <<"  is incosistent with node. \n.";
                    return SHIMMER_GENERIC_FAILURE;
                }

                //1.2.2 Sum inj
                for (size_t iC = 0; iC < NUM_GASES; iC++)
                {
                    double massfrac_face = massfrac_now(upw_num, iC);
                    rhs_nodes(iN,iC) += flux_inject * massfrac_face; 
                }
            }
        }
  

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

            if (node_prop.type == station_type::JUNCTION)
                continue;

            auto bnd =  node_prop.node_station->boundary();

            switch (node_prop.type) 
            {
                case station_type::ENTRY_L_REG: 
                case station_type::PRIVATE_INLET:
                {
                    auto v = std::abs(bnd.value(at_step));                
                    for (size_t iC = 0; iC < NUM_GASES; iC++)
                        rhs_nodes(node_prop.i_snum, iC) += v * massfrac_now(node_prop.i_snum, iC);
                    break;
                }
                case station_type::ENTRY_P_REG: 
                {
                    auto v = std::abs(var_all.L_rate(node_prop.i_snum));
                    for (size_t iC = 0; iC < NUM_GASES; iC++)
                        rhs_nodes(node_prop.i_snum, iC) += v * massfrac_now(node_prop.i_snum, iC);
                    break;
                }
                // Ejection
                case station_type::EXIT_L_REG: 
                case station_type::PRIVATE_OUTLET: 
                {
                    auto v = std::abs(bnd.value(at_step));                
                    lhs_nodes(node_prop.i_snum) += v;
                    break;        
                }
                default:
                    throw std::invalid_argument("QT Error: station is not valid.");
            }
        }
        
        if(unsteady)
        {
        // if NODE_ACCUMULATES
        // 1.4 Time term
        // V/c2 (dpdt) => I would rather do  V*(d(p/c2)/dt) thus I would need also c2 in t^{n} 

            vector_t phi = vector_t::Zero(infra_.num_original_stations);

            size_t count = 0; 
            for(auto itor = v_range.first; itor != v_range.second; itor++)
            {
                const auto & node_prop = infra_.graph[*itor]; 

                auto rho_now = rho_all_nodes_evol_(at_step,node_prop.i_snum);
                auto rho_old = rho_all_nodes_evol_(at_step-1,node_prop.i_snum);
                lhs_nodes(node_prop.i_snum) +=volume(*itor, infra_.graph) * (rho_now - rho_old) / dt;

                count++;
                if(count == infra_.num_original_stations)
                    break;
            }
        }
        // 3. 4 Solve Y^n+1            
        matrix_t lhs_inv =  lhs_nodes.cwiseInverse().asDiagonal();

        massfrac_next.topRows( num_nodes) = lhs_inv * rhs_nodes; 
        
        clip_and_renormalize(massfrac_next.topRows( num_nodes));
        
        return true;
    }


    void
    quality_tracking(double dt, 
                    const variable& var_all,
                    const matrix_t & massfrac_now, 
                    matrix_t& massfrac_next,
                    const vector_t & rho_all_previous, 
                    const vector_t & rho_nodes_previous)
    {
        /* 
            Solve: 

                dq/dt + d(v * q)/dx = 0 
            
            with 
                q = rho * Y      and      F = v * q
        */
        
        // Loop over network (original) pipes
        for(const auto& pd : infra_.pipe_discretizations)
        {
            auto dtdx =  dt/pd.dx; 
            auto num_loc_nodes = pd.nodelist.size();
            auto num_loc_pipes = pd.pipelist.size();

            ///1. Variables at fictitious nodes (primal mesh)/pipes (dual mesh) of the current pipe 
             
            matrix_t flux_nodes= massfracflux_fictitious_nodes(pd, var_all,
                                                              rho_all_previous,
                                                              massfrac_now);
            vector_t V_pipes   = velocity_fictitious_pipes(pd, var_all, 
                                                           rho_all_previous,
                                                           area_all_pipes_);
            vector_t rho_nodes = density_fictitious_nodes(pd,    
                                                          rho_nodes_previous);
            matrix_t q_nodes = matrix_t::Zero(num_loc_nodes, NUM_GASES); 

            for (int iN = 0; iN < num_loc_nodes; iN++)
            {
                auto idx = pd.nodelist[iN];
                q_nodes.row(iN) = rho_nodes[iN] * massfrac_now.row(idx);
            }
            
            // 2. Solve dq/dt + d(vq)/dx = 0
            for (int iN = 1; iN < num_loc_nodes-1; iN++)
            {
                auto idx = pd.nodelist[iN];

                // Values on interfaces of the primal mesh (ficts nodes) = dual mesh (fict pipes)
                auto [QL, QR] = LaxWendroff( dtdx, q_nodes, flux_nodes, iN);

                //auto [QF, QR] = MUSCL(); 

                vector_row_t FL = flux_exact(QL,V_pipes[iN-1]); 
                vector_row_t FR = flux_exact(QR,V_pipes[iN]);

                vector_row_t add = dtdx * (FR -FL); 
                vector_row_t q_next_i = q_nodes.row(iN) - dtdx * (FR -FL);

                // WARNING: massfrac_next = Qnext / rho_next...but I dont have rho_next 
                massfrac_next.row(idx) = q_next_i / rho_nodes[iN]; 

                clip_and_renormalize(massfrac_next.row(idx));

                #if 0 
                auto sum = 1.0 - massfrac_next.block(idx,0,1, NUM_GASES-1).sum(); 

                if(sum < -1.e-2)
                {
                    std::cout << "Problem in QT at node "<< iN << ": 1 - sum_{1}^{N-1} = " << sum << std::endl;
                    std::cout << massfrac_next << std::endl;
                    exit(1);
                }
                else if( std::abs(sum) < 1.e-2 )
                    sum = 0.;    

                massfrac_next(idx, NUM_GASES - 1) = sum;
                #endif

            }


        }
    }


    void
    fluid_dynamics(size_t it, 
                 double tol,
                 double dt,         // Most of this should be passed in a config struct 
                 size_t MAX_CONSTRAINT_ITER,
                 equation_of_state & eos)
    {
        bool unsteady = true;

        pipe_stations_activation(it, var_all_);

        size_t ic;
        for(ic = 0; ic <= MAX_CONSTRAINT_ITER; ic++)
        {
            std::cout<<"****************************************************************"<< std::endl;
            std::cout << " Iteration CONSTRAINTS it ..............."<<ic<< " ... at time "<< it << std::endl;
            std::cout<<"****************************************************************"<< std::endl;

            // To be finished when it is clear how molfrac changes and modifies mu.
            auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

            linearized_fluid_solver lfs(it, unsteady, tol, dt, temperature_, mu, inc_all_, infra_.graph);
            lfs.run(area_all_pipes_, var_all_guess_, var_all_, &eos, ic);

            bool pass_constr = lfs.check_constraints(it); 
            bool pass_control= lfs.check_controls(it);

            if(pass_constr && pass_control)
            {
                std::cout<< "++++++++++++++++++**** MODIFIED VARIABLE ****++++++++++++++++++++++ " << std::endl;
                var_all_ =  lfs.get_variable();
                rho_all_pipes_ =  eos.density(&lfs);  // needed for the computation of the velocity
                rho_all_nodes_ =  eos.density_nodes(&lfs);  // needed for the mass balance on network nodes 

                break;
            }
        }

        if(ic == MAX_CONSTRAINT_ITER)
            std::cerr << "ERROR: FAILURE to apply HARD constraints. Max number of iterations has been reached.";

        var_all_evol_.row(it) =  var_all_.make_vector();
        rho_all_pipes_evol_.row(it) =  rho_all_pipes_;
        rho_all_nodes_evol_.row(it) =  rho_all_nodes_;

        var_all_guess_ = var_all_;
    }


    bool
    steady_state(const variable& var_guess,
                 double dt,
                 double tolerance)
    {
        bool unsteady = false;
        bool converged = false; 

        EQ_OF_STATE eos;    
        auto [massfrac_nodes, massfrac_pipes] =  eos.molfrac_2_massfrac(infra_.graph, inc_all_);
        // Mass fraction by comp, needs total molar mass. So molar mass has to be updated any time molfrac changes!

        // Viscosity changes with molfrac (stored by comp in the graph).
        auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

        size_t iter = 0;
        auto var_time = variable(vector_t::Zero(num_vertices(infra_.graph)),
                                 vector_t::Zero(num_edges(infra_.graph)),
                                 vector_t::Zero(num_vertices(infra_.graph)));

        pipe_stations_activation(iter, var_time);

        variable var_previous = var_guess;
        variable var_current;

        double t = 0;
        for(size_t it = 0; it < MAX_ITERS_STEADY_; it++)
        {
            // 1. Update molar masses (mm) inside eos and molar frac inside graph
            auto [molfrac_nodes, molfrac_pipes] =  eos.massfrac_2_molfrac(massfrac_nodes, massfrac_pipes);

            auto index = get(boost::vertex_index, infra_.graph);
            for (auto vp = vertices(infra_.graph); vp.first != vp.second; ++vp.first) {
                auto idx = index[*vp.first]; 
                infra_.graph[idx].gas_mixture = molfrac_nodes.row(idx).transpose();
            }

            linearized_fluid_solver lfs(iter, unsteady, tolerance, dt, temperature_, mu, inc_all_, infra_.graph);
            lfs.run(area_all_pipes_, var_previous, var_time, &eos);
            var_current = lfs.get_variable();

            // 3. 1. Continuity at network nodes and fictitious nodes
            matrix_t massfrac_next_nodes = matrix_t::Zero(boost::num_vertices(infra_.graph) , NUM_GASES);

            admixing(unsteady, it, dt, 
                         var_current, 
                         massfrac_nodes, massfrac_next_nodes);


            auto error_pressure = (var_previous.pressure - var_current.pressure).norm(); 
            auto error_flux     = (var_previous.flux - var_current.flux).norm(); 
            auto error_massfrac = (massfrac_nodes - massfrac_next_nodes).norm();

            auto norm_pressure = (var_previous.pressure).norm(); 
            auto norm_flux   = (var_current.flux).norm(); 
            auto norm_massfrac = massfrac_nodes.norm();
            auto residual =  std::max(std::max(error_flux, error_pressure), error_massfrac) 
                               / std::max(std::max(norm_flux, norm_pressure), norm_massfrac); 

            if(residual < tolerance)
            {
                var_all_evol_.row(0) =  lfs.get_variable().make_vector();
                rho_all_pipes_evol_.row(0) =  eos.density(&lfs);
                rho_all_nodes_evol_.row(0) =  eos.density_nodes(&lfs);

                massfrac_pipes = inc_all_.matrix_in().transpose() * massfrac_nodes;    
                auto pair_molfrac = eos.massfrac_2_molfrac(massfrac_nodes, massfrac_pipes);
                molfrac_all_evol_[0] = pair_molfrac.first; 

                return true;
            }
            massfrac_nodes = massfrac_next_nodes;
            massfrac_pipes = inc_all_.matrix_in().transpose() * massfrac_nodes;    
            var_previous = var_current;
        }

        std::cerr << "Admixing - Steady has NOT CONVERGED." << std::endl;
        return SHIMMER_GENERIC_FAILURE;
    }

    void
    advance(double dt,        
            double tol)
    {
        bool unsteady = true;
        size_t MAX_CONSTRAINT_ITER = 10;

        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 

        matrix_t massfrac_now_nodes = massfrac_guess_;
        matrix_t massfrac_next_nodes = matrix_t::Zero(num_nodes, NUM_GASES);

        std::cout << std::setprecision(4 )<<" * massfrac_now : \n" << massfrac_guess_ << std::endl;

        molfrac_all_evol_[0] = build_molfrac_nodes(infra_.graph);

        EQ_OF_STATE eos;

        double t = 0;
        for(size_t it = 1; it < num_time_steps_; it++, t+=dt)
        {
            std::cout<<"========================================================"<< std::endl;
            std::cout << "Solving at time ...."<< it <<std::endl;
            std::cout<<"========================================================"<< std::endl;

            // 1. Update molar masses (mm) inside eos and molar frac inside graph
            matrix_t massfrac_now_pipes = inc_all_.matrix_in().transpose() * massfrac_now_nodes;    
            auto [molfrac_nodes, molfrac_pipes] =  eos.massfrac_2_molfrac(massfrac_now_nodes, massfrac_now_pipes);

            auto index = get(boost::vertex_index, infra_.graph);
            for (auto vp = vertices(infra_.graph); vp.first != vp.second; ++vp.first) {
                auto idx = index[*vp.first]; 
                infra_.graph[idx].gas_mixture = molfrac_nodes.row(idx).transpose();
            }

            // 2. Fluid solver
            fluid_dynamics(it, tol, dt, MAX_CONSTRAINT_ITER, eos);

            // 3. [massfrac_nodes, massfrac_pipes] = quality_tracking();
            //matrix_t lhs_nodes = vector_t::Zero(infra_.num_original_stations, NUM_GASES); 
            vector_t lhs_nodes = vector_t::Zero(infra_.num_original_stations); 
            matrix_t rhs_nodes = matrix_t::Zero(infra_.num_original_pipes, NUM_GASES); 
            massfrac_next_nodes = matrix_t::Zero(num_nodes, NUM_GASES);

            // 3. 1. Continuity at network nodes
            admixing(unsteady, it, dt, 
                         var_all_, 
                         massfrac_now_nodes, massfrac_next_nodes);

            if(refine_)
            {
                // 3. 2. Continuity at fictitious nodes 
                quality_tracking(dt, var_all_, 
                                 massfrac_now_nodes, 
                                 massfrac_next_nodes,
                                 rho_all_pipes_evol_.row(it-1),
                                 rho_all_nodes_evol_.row(it-1) );
            }

            massfrac_now_nodes = massfrac_next_nodes;
            molfrac_all_evol_[it] = build_molfrac_nodes(infra_.graph); 
            //std::cout << std::setprecision(4 )<<" * massfrac_now : \n" << massfrac_now_nodes << std::endl;
        }
        return;
    }

/*
    std::pair<vector_t, vector_t>
    MUSCL(double dtdx, const matrix_t& Q, size_t i)
    {
        auto kappa = 1.0/3.0;

        vector_t QL, QR;

        if(i > 2)        
           QL =  Q.row(i-1) + 0.25 * (1.0 -kappa) * (Q.row(i-1) + Q.row(i-2)) +  0.25 * (1.0 +kappa) * (Q.row(i) + Q.row(i-1));
        else
           QL =  Q.row(i-1) + 0.25 * (1.0 -kappa) * (Q.row(i-1) + Q.row(i-2)) +  0.25 * (1.0 +kappa) * (Q.row(i) + Q.row(i-1));

        if()           
           QR =  Q.row(i) + 0.25 * (1.0 -kappa) * (Q.row(i) + Q.row(i-1)) +  0.25 * (1.0 +kappa) * (Q.row(i+1) + Q.row(i));

        return std::make_pair(QL, QR);
    }
*/ 

    matrix_t
    massfracflux_fictitious_nodes(const pipe_discretization& pd, 
                                  const variable& var_all, 
                                  const vector_t& rho_all, 
                                  const matrix_t& massfrac_now)
    {
        auto num_loc_nodes = pd.nodelist.size();
        auto num_loc_pipes = pd.pipelist.size();

        // Mass flux within fictitious pipes(dual mesh):
        //      
        //       F_mass_j = rho_j * vel_j
        // 
        vector_t Fmass = massflux_fictitious_pipes(pd, var_all, area_all_pipes_);

        assert(Fmass.size()+1  == num_loc_nodes && "Incorrect size for local velocities");

        // 1. Fluxes on fictitius nodes (primal mesh)
        //      
        //       F_i = rho_i * vel_i * Y_i 
        // 
        matrix_t F = matrix_t::Zero(num_loc_nodes, NUM_GASES); 
        for (int iN = 1; iN < num_loc_nodes-1; iN++)
        {
            auto idx = pd.nodelist[iN];
            F.row(iN) = 0.5 * (Fmass(iN) + Fmass(iN-1)) * massfrac_now.row(idx);
        }
        // outer points: chosen to be equal to the network volume value. 
        F.row(0) =  Fmass(0) * massfrac_now.row(pd.nodelist[0]);
        F.row(num_loc_nodes -1) =  Fmass(num_loc_pipes-1) * massfrac_now.row(pd.nodelist[num_loc_nodes -1]);

        return F;
    }

    vector_t 
    density_fictitious_nodes(const pipe_discretization& pd, const vector_t& rho_nodes_all)
    {
        vector_t rho =  vector_t::Zero(pd.nodelist.size());

        for (int iN = 0; iN < pd.nodelist.size(); iN++)
        {
            auto idx = pd.nodelist[iN];
            rho(iN) = rho_nodes_all(idx);
        }
        return rho;
    }

    vector_t
    massflux_fictitious_pipes(const pipe_discretization& pd, const variable& var, const vector_t& area) const
    {
        vector_t mf = vector_t::Zero(pd.pipelist.size());
         
        for (auto iCount = 0;  iCount < pd.pipelist.size(); iCount++)
        {
            auto pipe_num = pd.pipelist[iCount];
            mf(iCount) = var.flux(pipe_num) / area(pipe_num);
        }

        return mf;      
    }

    vector_t
    velocity_fictitious_pipes(const pipe_discretization& pd,
                              const variable& var,
                              const vector_t& rho_pipes,
                              const vector_t& area_pipes) const
    {
        vector_t vel = vector_t::Zero(pd.pipelist.size());
         
        for (auto iCount = 0;  iCount < pd.pipelist.size(); iCount++)
        {
            auto pipe_num = pd.pipelist[iCount];
            vel(iCount) = var.flux(pipe_num) 
                                 / (rho_pipes(pipe_num) * area_pipes(pipe_num));
        }

        return vel;      
    }

    matrix_t 
    velocity_evolution() const
    {
        // This one cannot be called when doing only steady state 

        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 
        /// vel [m/s] velocity of the gas within pipes.
        matrix_t flux_evol = var_all_evol_.middleCols(num_nodes, num_pipes);
        matrix_t arho_evol = rho_all_pipes_evol_;

        for(auto iP = 0; iP < num_pipes; iP++)
            arho_evol.col(iP) *= area_all_pipes_(iP);  
        return flux_evol.array() / arho_evol.array();       
    }

    vector_t solution() const {return var_all_.make_vector();}
    vector_t guess() const {return var_all_guess_.make_vector();}
    matrix_t solution_evolution() const{ return var_all_evol_; }
    std::vector<matrix_t> molar_fractions_evolution() const{ return molfrac_all_evol_;}

};

}
