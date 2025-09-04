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

#include "solver/gas_law.h"
#include <iomanip>

namespace shimmer{


matrix_t
build_x_nodes(const infrastructure_graph& g)
{
    size_t num_nodes = boost::num_vertices(g);
    Eigen::MatrixXd x(num_nodes, NUM_GASES);

    auto index = get(boost::vertex_index, g);

    for (auto vp = vertices(g); vp.first != vp.second; ++vp.first) {
        auto idx = index[*vp.first]; 
        x.row(idx) = g[idx].gas_mixture;
    }

    return x;
}


vector_t 
equation_of_state::density(linearized_fluid_solver *lfs) const
{
    auto [c2_nodes, c2_pipes] = speed_of_sound(lfs);

    return lfs->pressure_pipes().array() / c2_pipes.array();
}


vector_t 
equation_of_state::density_nodes(linearized_fluid_solver *lfs) const
{
    auto [c2_nodes, c2_pipes] = speed_of_sound(lfs);

    return lfs->pressure_nodes().array() / c2_nodes.array();
}


vector_t
equation_of_state::compute_R(const vector_t& molar_mass) const
{
    return Runiversal_ * molar_mass.cwiseInverse();     
}

void
equation_of_state::mixture_molar_mass(const infrastructure_graph& graph, const incidence& inc)
{
    // Molar frac by comp and by pipe/node
    matrix_t x_nodes = build_x_nodes(graph);
    matrix_t x_pipes = inc.matrix_in().transpose() * x_nodes;   

    mixture_molar_mass(x_nodes, x_pipes);
}


// ----------------------------------------------------------------------------
//                    Papay equation of state
// ----------------------------------------------------------------------------

papay::papay(): T_cr_(190.6), p_cr_(4.595e+6)
{}


void 
papay::initialization(linearized_fluid_solver *lfs){}


void
papay::mixture_molar_mass(const matrix_t& x_nodes, const matrix_t& x_pipes)
{
    mm_nodes_ = vector_t::Zero(x_nodes.rows()); 
    mm_pipes_ = vector_t::Zero(x_pipes.rows()); 

    mm_nodes_.setConstant(16.0);
    mm_pipes_.setConstant(16.0);
}


vector_t 
papay::compute_Z(double temperature, const vector_t& pressure) const
{
    auto coef0 = std::exp(-2.260 * (temperature/T_cr_)) * (3.52 / p_cr_);
    auto coef1 = std::exp(-1.878 * (temperature/T_cr_)) * 0.274/(p_cr_*p_cr_); 
    return 1.0 - coef0 *  pressure.array() 
                        + coef1 * pressure.array().square() ;

    //return Z;
}


std::pair<vector_t, vector_t>
papay::speed_of_sound(linearized_fluid_solver *lfs) const
{
    vector_t Z_nodes = compute_Z(lfs->temperature(), lfs->pressure_nodes());
    vector_t Z_pipes = compute_Z(lfs->temperature(), lfs->pressure_pipes());

    vector_t R_nodes = compute_R(mm_nodes_); 
    vector_t R_pipes = compute_R(mm_pipes_); 

    vector_t c2_nodes = Z_nodes.cwiseProduct(R_nodes) * lfs->temperature(); 
    vector_t c2_pipes = Z_pipes.cwiseProduct(R_pipes) * lfs->temperature();

    return std::make_pair(c2_nodes, c2_pipes);
}


std::pair<matrix_t, matrix_t>
papay::molarfrac_2_massfrac(const infrastructure_graph& graph, const incidence& inc)
{
    matrix_t y_nodes = matrix_t::Zero(boost::num_vertices(graph), NUM_GASES ); 
    matrix_t y_pipes = matrix_t::Zero(boost::num_edges(graph), NUM_GASES ); 

    y_nodes.col(0) = vector_t::Ones(boost::num_vertices(graph), 1); 
    y_pipes.col(0) = vector_t::Ones(boost::num_edges(graph), 1 ); 

    //In order to not define more variable, we use y as input, instead of x, since they are equal. 
    mixture_molar_mass(y_nodes, y_pipes);

    return std::make_pair(y_nodes, y_pipes);
}


std::pair<matrix_t, matrix_t> 
papay::massfrac_2_molarfrac(const matrix_t& y_nodes, const matrix_t& y_pipes)
{
    matrix_t x_nodes = matrix_t::Zero(y_nodes.rows(), NUM_GASES ); 
    matrix_t x_pipes = matrix_t::Zero(y_pipes.rows(), NUM_GASES ); 

    x_nodes.col(0) = vector_t::Ones(y_nodes.rows());
    x_pipes.col(0) = vector_t::Ones(y_pipes.rows() ); 

    mixture_molar_mass(x_nodes, x_pipes);
    
    return std::make_pair(y_nodes, y_pipes);
}
// ----------------------------------------------------------------------------
//                    GERG equation of state
// ----------------------------------------------------------------------------

#ifdef HAVE_MATLAB_GERG
gerg_params::gerg_params(const gerg_reducing_params_t&  rp, 
                const gerg_pseudo_critical_pt_t& psc,
                const gerg_thermo_params_t& tp): 
                reducing_params(rp),
                pseudo_critical_pt (psc),
                params (tp) {}; 


gerg_params:: gerg_params(){};


gerg::gerg()
{
    mm_component = vector_t::Zero(21); 

    mm_component <<
    //  CH4              N2             CO2             C2H6
    1.604246e+01,   2.801340e+01,   4.400950e+01,   3.006904e+01,
    //  C3H8          i-C4H10         n-C4H10         i-C5H12
    4.409562e+01,   5.812220e+01,   5.812220e+01,   7.214878e+01,
    //n-C5H12           C6H14           C7H16           C8H18
    7.214878e+01,   8.617536e+01,   1.0020194e+02,  1.1422852e+02,
    //  C9H20           C10H22          H2              O2
    1.282551e+02,   1.4228168e+02,  2.015880e+00,   3.199880e+01,
    //  CO              H2O             H2S             He
    2.801010e+01,   1.801528e+01,   3.408088e+01,   4.002602e+00,
    //  Ar
    3.9948000e+01;

    gas_name_ = std::vector<int>
       {GAS_TYPE::CH4, GAS_TYPE::N2, GAS_TYPE::CO2, GAS_TYPE::C2H6,
        GAS_TYPE::C3H8, GAS_TYPE::i_C4H10, GAS_TYPE::n_C4H10,
        GAS_TYPE::i_C5H12, GAS_TYPE::n_C5H12,GAS_TYPE::C6H14,
        GAS_TYPE::C7H16, GAS_TYPE::C8H18, GAS_TYPE::C9H20,
        GAS_TYPE::C10H22,GAS_TYPE::H2, GAS_TYPE::O2, GAS_TYPE::CO,
        GAS_TYPE::H2O,GAS_TYPE::H2S,GAS_TYPE::He,GAS_TYPE::Ar}; 

}


gerg_params
gerg::compute_params(const matrix_t& x)
{
    auto rp  = GERG::reducing_parameters(x);
    auto psc = GERG::pseudo_critical_point(x, x.rows());
    auto tp  = gerg_thermo_params_t();
    tp.Type = gerg_thermo_params_t::Types::Gas_phase;

    return gerg_params(rp, psc, tp); 

}


void
gerg::initialization(linearized_fluid_solver *lfs)
{
    params_nodes_ = compute_params(lfs->x_nodes());
    params_pipes_ = compute_params(lfs->x_pipes());
}


vector_t 
gerg::compute_Z( const double& temperature,
                 const vector_t& pressure,
                 const matrix_t& x,
                 const gerg_params& gp) const
{
    auto eos =  GERG::thermodynamic_properties(pressure, temperature, x,
                                              gp.reducing_params,
                                              gp.pseudo_critical_pt,
                                              gp.params);
    return eos.Z;
}


std::pair<vector_t, vector_t>
gerg::speed_of_sound(linearized_fluid_solver *lfs) const
{
    vector_t Z_nodes = compute_Z(lfs->temperature(), lfs->pressure_nodes(), lfs->x_nodes(), params_nodes_);
    vector_t Z_pipes = compute_Z(lfs->temperature(), lfs->pressure_pipes(), lfs->x_pipes(), params_pipes_);

    vector_t R_nodes = compute_R(mm_nodes_); 
    vector_t R_pipes = compute_R(mm_pipes_); 

    vector_t c2_nodes = Z_nodes.cwiseProduct(R_nodes) * lfs->temperature(); 
    vector_t c2_pipes = Z_pipes.cwiseProduct(R_pipes) * lfs->temperature();

    return std::make_pair(c2_nodes, c2_pipes);
}


void
gerg::mixture_molar_mass(const matrix_t& x_nodes, const matrix_t& x_pipes)
{
    mm_nodes_ = vector_t::Zero(x_nodes.rows()); 
    mm_pipes_ = vector_t::Zero(x_pipes.rows()); 

    for(size_t i = 0; i <= 20; i++)
        mm_nodes_ +=  mm_component(i) * x_nodes.col(gas_name_[i]); 

    for(size_t i = 0; i <= 20; i++)
        mm_pipes_ +=  mm_component(i) * x_pipes.col(gas_name_[i]); 
}


std::pair<matrix_t, matrix_t>
gerg::molarfrac_2_massfrac(const infrastructure_graph& graph, const incidence& inc)
{
    // Molar frac by comp and by pipe/node
    matrix_t x_nodes = build_x_nodes(graph);
    matrix_t x_pipes = inc.matrix_in().transpose() * x_nodes;      
    mixture_molar_mass(x_nodes, x_pipes);

    matrix_t y_nodes = matrix_t::Zero(boost::num_vertices(graph), NUM_GASES ); 
    matrix_t y_pipes = matrix_t::Zero(boost::num_edges(graph), NUM_GASES ); 

    for(size_t iN = 0; iN < boost::num_vertices(graph); iN++)
        y_nodes.row(iN) =  mm_component.array() * x_nodes.row(iN).array() / mm_nodes_.array(); 

    for(size_t iP = 0; iP < boost::num_edges(graph); iP++)
        y_pipes.row(iP) =  mm_component.array() * x_pipes.row(iP).array() / mm_pipes_.array(); 

    return std::make_pair(y_nodes, y_pipes);
}

std::pair<matrix_t, matrix_t> 
gerg::massfrac_2_molarfrac(const matrix_t&y_nodes, const matrix_t&y_pipes)
{
    matrix_t x_nodes = matrix_t::Zero(y_nodes.rows(), NUM_GASES ); 
    matrix_t x_pipes = matrix_t::Zero(y_pipes.rows(), NUM_GASES ); 

    for(int iN = 0; iN < NUM_GASES; iN++)
    {
        //diag
        matrix_t A =  mm_component.asDiagonal();

        // off diag
        for(int iComp = 0; iComp < NUM_GASES; iComp++)
        {
            A.row(iComp) -= y_nodes(iN, iComp) * mm_component;
        }

        x_nodes.row(iN) = A.fullPivLu().solve(vector_t::Zero(NUM_GASES)).transpose();
    }
    
    // Update mixture molar mass
    mixture_molar_mass(x_nodes, x_pipes);

    return std::make_pair(x_nodes, x_pipes);
}

#endif /* HAVE_MATLAB_GERG */
// ----------------------------------------------------------------------------
//                    GERG equation of state AGA8CODE
// ----------------------------------------------------------------------------

gerg_aga::gerg_aga()
{
    tolerance_ = 1.e-12;

    shimmer_gerg::gerg_functions::setup_GERG();

    mm_component = shimmer_gerg::gerg_functions::component_molar_masses(); 

}


void
gerg_aga::initialization(linearized_fluid_solver *lfs)
{}


vector_t 
gerg_aga::compute_Z( const vector_t& temperature,
                   const vector_t& pressure,
                   const matrix_t& x) const
{
    auto type = gerg_aga_thermo_params_t::Types::Gas_phase;
    
    auto eos =  shimmer_gerg::gerg_functions::thermodynamic_properties(
                        temperature,
                        pressure, //[Pa]
                        x,
                        type,
                        tolerance_);
    return eos.Z;
}


std::pair<vector_t, vector_t>
gerg_aga::speed_of_sound(linearized_fluid_solver *lfs) const
{
    vector_t T_nodes = vector_t::Zero(lfs->num_nodes());
    vector_t T_pipes = vector_t::Zero(lfs->num_pipes());

    T_nodes.setConstant(lfs->temperature());
    T_pipes.setConstant(lfs->temperature());

    vector_t Z_nodes = compute_Z( T_nodes, lfs->pressure_nodes(), lfs->x_nodes());
    vector_t Z_pipes = compute_Z( T_pipes, lfs->pressure_pipes(), lfs->x_pipes());
    vector_t R_nodes = compute_R(mm_nodes_); 
    vector_t R_pipes = compute_R(mm_pipes_); 

    vector_t c2_nodes = Z_nodes.array() * R_nodes.array() * lfs->temperature(); 
    vector_t c2_pipes = Z_pipes.array() * R_pipes.array() * lfs->temperature();

    return std::make_pair(c2_nodes, c2_pipes);
}


void
gerg_aga::mixture_molar_mass(const matrix_t& x_nodes, const matrix_t& x_pipes)
{
    mm_nodes_ = vector_t::Zero(x_nodes.rows()); 
    mm_pipes_ = vector_t::Zero(x_pipes.rows()); 

    shimmer_gerg::gerg_functions::molar_mass(x_nodes, tolerance_, mm_nodes_);
    shimmer_gerg::gerg_functions::molar_mass(x_pipes, tolerance_, mm_pipes_);
}


std::pair<matrix_t, matrix_t>
gerg_aga::molarfrac_2_massfrac(const infrastructure_graph& graph, const incidence& inc)
{
    // Molar frac by comp and by pipe/node
    // Molar frac by comp and by pipe/node
    matrix_t x_nodes = build_x_nodes(graph);
    matrix_t x_pipes = inc.matrix_in().transpose() * x_nodes;      
    mixture_molar_mass(x_nodes, x_pipes);

    matrix_t y_nodes = matrix_t::Zero(boost::num_vertices(graph), NUM_GASES ); 
    matrix_t y_pipes = matrix_t::Zero(boost::num_edges(graph), NUM_GASES ); 
      
    for(size_t iN = 0; iN < boost::num_vertices(graph); iN++)
        y_nodes.row(iN) =  mm_component.transpose().array() * x_nodes.row(iN).array() / mm_nodes_(iN); 

    for(size_t iP = 0; iP < boost::num_edges(graph); iP++)
        y_pipes.row(iP) =  mm_component.transpose().array() * x_pipes.row(iP).array() / mm_pipes_(iP); 

    return std::make_pair(y_nodes, y_pipes);
}

std::pair<matrix_t, matrix_t> 
gerg_aga::massfrac_2_molarfrac(const matrix_t&y_nodes, const matrix_t&y_pipes)
{
    matrix_t x_nodes = matrix_t::Zero(y_nodes.rows(), NUM_GASES ); 
    matrix_t x_pipes = matrix_t::Zero(y_pipes.rows(), NUM_GASES ); 

    for(int iN = 0; iN < y_nodes.rows(); iN++)
    {
        //diag
        matrix_t A =  mm_component.asDiagonal();

        // off diag
        for(int iComp = 0; iComp < NUM_GASES; iComp++)
        {
            A.row(iComp) -= y_nodes(iN, iComp) * mm_component;
        }

        x_nodes.row(iN) = A.fullPivLu().solve(vector_t::Zero(NUM_GASES)).transpose();
    }
    
    // Update molar mixture
    mixture_molar_mass(x_nodes, x_pipes);

    return std::make_pair(x_nodes, x_pipes);
}

} //end namespace shimmer
