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
    int MAX_ITERS_STEADY_;

    matrix_t massfrac_guess_;

    vector_t rho_msh_;
    vector_t rho_nodes_;
    variable var_msh_;
    variable var_msh_guess_;
    vector_t area_msh_pipes_;

    matrix_t rho_msh_in_time_;
    matrix_t rho_nodes_in_time_;
    matrix_t var_msh_in_time_;
    std::vector<matrix_t> molfrac_in_time_;

    incidence inc_msh_;
    infrastructure& infra_;

public:
    qt_solver(infrastructure& infrain,
              double Tm,
              const bool& refine):
              infra_(infrain), temperature_(Tm), refine_(refine)
    {
        MAX_ITERS_STEADY_ = 500;
        inc_msh_ = incidence(infra_.graph);
        area_msh_pipes_ = area(infra_.graph);
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
        rho_nodes_ = vector_t::Ones(boost::num_vertices(infra_.graph));
    }


    void
    initialization( const variable& var_guess,
                    double dt,
                    double tolerance)
    {
        bool unsteady = false;

        EQ_OF_STATE eos;    
        auto [massfrac_nodes_, massfrac_pipes_] =  eos.molfrac_2_massfrac(infra_.graph, inc_msh_);
        // Mass fraction by comp, needs total molar mass. So molar mass has to be updated any time molfrac changes!

        // Viscosity changes with molfrac (stored by comp in the graph).
        auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

        size_t iter = 0;
        auto var_time = variable(vector_t::Zero(num_vertices(infra_.graph)),
                                 vector_t::Zero(num_edges(infra_.graph)),
                                 vector_t::Zero(num_vertices(infra_.graph)));

        pipe_stations_activation(iter, var_time);

        linearized_fluid_solver lfs(iter, unsteady, tolerance, dt, temperature_, mu, inc_msh_, infra_.graph);
        lfs.run(area_msh_pipes_, var_guess, var_time, &eos);
        var_msh_guess_ = lfs.get_variable();
        massfrac_guess_ = massfrac_nodes_;
        rho_msh_ = eos.density(&lfs);
        rho_nodes_ = eos.density_nodes(&lfs);
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
    qt_network_nodes(bool unsteady, size_t at_step, double dt,  const variable& var_msh,
                     const  matrix_t& massfrac_now, matrix_t& massfrac_next)
    {
        vector_t lhs_nodes = vector_t::Zero(infra_.num_original_stations);
        matrix_t rhs_nodes = matrix_t::Zero(infra_.num_original_stations, NUM_GASES );

        // Mass conservation for network nodes
        auto inc_mat = inc_msh_.matrix();


        // 1. Pipes Injection/Ejection
        
        using svec_itor_t = Eigen::SparseVector<double>::InnerIterator;

        // Here, local pipe quantities are needed for network(original) nodes
        auto inc_smat = inc_msh_.matrix();
        for(size_t iN = 0; iN < infra_.num_original_stations; iN++)
        {
            Eigen::SparseVector<double> inc_node_i = inc_msh_.matrix().row(iN);
            
            // Loop by face
            for(svec_itor_t it(inc_node_i); it; ++it)
            {
                auto  pipe_num = it.row();

                double inc_val = it.value();
                double val = inc_val * var_msh.flux(pipe_num); 

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
                    auto v = std::abs(var_msh.L_rate(node_prop.i_snum));
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

                auto rho_now = rho_nodes_in_time_(at_step,node_prop.i_snum);
                auto rho_old = rho_nodes_in_time_(at_step-1,node_prop.i_snum);
                lhs_nodes(node_prop.i_snum) +=volume(*itor, infra_.graph) * (rho_now - rho_old) / dt;

                count++;
                if(count == infra_.num_original_stations)
                    break;
            }
        }
        // 3. 4 Solve Y^n+1            
        matrix_t lhs_inv =  lhs_nodes.cwiseInverse().asDiagonal();

        massfrac_next.topRows( infra_.num_original_stations) = lhs_inv * rhs_nodes; 

        return true;
    }


    void
    qt_fictitious_nodes(double dt, const variable& var_msh,
                    const matrix_t& massfrac_now, matrix_t& massfrac_next,
                    const vector_t& rho_msh_previous,const vector_t & rho_nodes_previous)
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
             
            matrix_t flux_nodes= massfracflux_fictitious_nodes(pd, var_msh,
                                                              rho_msh_previous,
                                                              massfrac_now);
            vector_t V_pipes   = velocity_fictitious_pipes(pd, var_msh, 
                                                           rho_msh_previous,
                                                           area_msh_pipes_);
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
            }


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

            // To be finished when it is clear how molfrac changes and modifies mu.
            auto mu = viscosity<viscosity_type>(temperature_, infra_.graph);

            linearized_fluid_solver lfs(it, unsteady, tol, dt, temperature_, mu, inc_msh_, infra_.graph);
            lfs.run(area_msh_pipes_, var_msh_guess_, var_msh_, &eos, ic);

            bool pass_constr = lfs.check_constraints(it); 
            bool pass_control= lfs.check_controls(it);

            if(pass_constr && pass_control)
            {
                std::cout<< "++++++++++++++++++**** MODIFIED VARIABLE ****++++++++++++++++++++++ " << std::endl;
                var_msh_ =  lfs.get_variable();
                rho_msh_ =  eos.density(&lfs);  // needed for the computation of the velocity
                rho_nodes_ =  eos.density_nodes(&lfs);  // needed for the mass balance on network nodes 

                break;
            }
        }

        if(ic == MAX_CONSTRAINT_ITER)
            std::cout << "ERROR: FAILURE to apply HARD constraints. Max number of iterations has been reached.";

        var_msh_in_time_.row(it) =  var_msh_.make_vector();
        rho_msh_in_time_.row(it) =  rho_msh_;
        rho_nodes_in_time_.row(it) =  rho_nodes_;

        var_msh_guess_ = var_msh_;
    }


    void
    advance(double dt,
            size_t num_steps,            
            double tol)
    {
        bool unsteady = true;
        size_t MAX_CONSTRAINT_ITER = 10;

        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 

        var_msh_in_time_ = matrix_t::Zero(num_steps, num_pipes + 2 * num_nodes);
        rho_msh_in_time_ = matrix_t::Ones(num_steps, num_pipes);
        rho_nodes_in_time_ = matrix_t::Ones(num_steps, num_nodes);
        molfrac_in_time_ = std::vector<matrix_t>(num_steps);

        var_msh_ = var_msh_guess_;
        var_msh_in_time_.row(0) =  var_msh_.make_vector();
        rho_msh_in_time_.row(0) =  rho_msh_.transpose();
        rho_nodes_in_time_.row(0) =  rho_nodes_.transpose();

        matrix_t massfrac_now_nodes = massfrac_guess_;
        matrix_t massfrac_next_nodes = matrix_t::Zero(num_nodes, NUM_GASES);

        std::cout << std::setprecision(4 )<<" * massfrac_now :" << massfrac_guess_ << std::endl;

        molfrac_in_time_[0] = build_molfrac_nodes(infra_.graph);

        EQ_OF_STATE eos;

        double t = 0;
        for(size_t it = 1; it < num_steps; it++, t+=dt)
        {
            std::cout<<"========================================================"<< std::endl;
            std::cout << "Solving at time ...."<< it <<std::endl;
            std::cout<<"========================================================"<< std::endl;

            std::cout << std::setprecision(16 ) <<  "rho = [" << rho_nodes_in_time_.row(it-1) << "]; "<<std::endl;

            // 1. Update molar masses (mm) inside eos and molar frac inside graph
            matrix_t massfrac_now_pipes = inc_msh_.matrix_in().transpose() * massfrac_now_nodes;    
            auto [molfrac_nodes, molfrac_pipes] =  eos.massfrac_2_molfrac(massfrac_now_nodes, massfrac_now_pipes);

            auto index = get(boost::vertex_index, infra_.graph);
            for (auto vp = vertices(infra_.graph); vp.first != vp.second; ++vp.first) {
                auto idx = index[*vp.first]; 
                infra_.graph[idx].gas_mixture = molfrac_nodes.row(idx).transpose();
            }

            // 2. Fluid solver
            transmission(it, tol, dt, MAX_CONSTRAINT_ITER, eos);

            // 3. [massfrac_nodes, massfrac_pipes] = quality_tracking();
            //matrix_t lhs_nodes = vector_t::Zero(infra_.num_original_stations, NUM_GASES); 
            vector_t lhs_nodes = vector_t::Zero(infra_.num_original_stations); 
            matrix_t rhs_nodes = matrix_t::Zero(infra_.num_original_pipes, NUM_GASES); 

            // 3. 1. Continuity at network nodes 
            qt_network_nodes(unsteady, it, dt, 
                         var_msh_, 
                         massfrac_now_nodes, massfrac_next_nodes);

            if(refine_)
            {
                // 3. 2. Continuity at fictitious nodes 
                qt_fictitious_nodes( dt, var_msh_, 
                                     massfrac_now_nodes, 
                                     massfrac_next_nodes,
                                     rho_msh_in_time_.row(it-1),
                                     rho_nodes_in_time_.row(it-1) );
            }

            massfrac_now_nodes = massfrac_next_nodes;
            
            std::cout << std::setprecision(4 )<<" * massfrac_now :" << massfrac_now_nodes << std::endl;

            molfrac_in_time_[it] = build_molfrac_nodes(infra_.graph); 

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
    massfracflux_fictitious_nodes(const pipe_discretization& pd, const variable& var_msh, const vector_t& rho_msh, 
                        const matrix_t& massfrac_now)
    {
        auto num_loc_nodes = pd.nodelist.size();
        auto num_loc_pipes = pd.pipelist.size();

        // Mass flux within fictitious pipes(dual mesh):
        //      
        //       F_mass_j = rho_j * vel_j
        // 
        vector_t Fmass = massflux_fictitious_pipes(pd, var_msh, area_msh_pipes_);

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
    velocity_full() const
    {
        // This one cannot be called when doing only steady state 

        auto num_nodes = num_vertices(infra_.graph); 
        auto num_pipes = num_edges(infra_.graph); 
        /// vel [m/s] velocity of the gas within pipes.
        matrix_t flux_in_time = var_msh_in_time_.middleCols(num_nodes, num_pipes);
        matrix_t arho_in_time = rho_msh_in_time_;

        for(auto iP = 0; iP < num_pipes; iP++)
            arho_in_time.col(iP) *= area_msh_pipes_(iP);  
        return flux_in_time.array() / arho_in_time.array();       
    }

    vector_t solution() const {return var_msh_.make_vector();}
    vector_t guess() const {return var_msh_guess_.make_vector();}
    matrix_t solution_full() const{ return var_msh_in_time_; }
    std::vector<matrix_t> molar_fractions_full() const{ return molfrac_in_time_;}

};

}
