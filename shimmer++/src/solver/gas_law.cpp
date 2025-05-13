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

void
equation_of_state::compute_density(linearized_fluid_solver *lfs,
                                   const vector_t& c2_pipes)
{
    density_ = lfs->pressure_pipes().array() / c2_pipes.array();
}


const vector_t& 
equation_of_state::density()
{
    return density_; 
}

// ----------------------------------------------------------------------------
//                    Papay equation of state
// ----------------------------------------------------------------------------

papay::papay(): T_cr_(190.6), p_cr_(4.595e+6)
{}


void 
papay::initialization(linearized_fluid_solver *lfs){}



void
papay::compute_molar_mass(const matrix_t& y_nodes, const matrix_t& y_pipes)
{
    mm_nodes_ = vector_t::Zero(y_nodes.rows()); 
    mm_pipes_ = vector_t::Zero(y_pipes.rows()); 

    mm_nodes_.setConstant(16.0);
    mm_pipes_.setConstant(16.0);

    R_nodes_= Runiversal_ * mm_nodes_.cwiseInverse();     
    R_pipes_= Runiversal_ * mm_pipes_.cwiseInverse();     
}

vector_t 
papay::compute(double temperature, const vector_t& pressure)
{
    auto coef0 = std::exp(-2.260 * (temperature/T_cr_)) * (3.52 / p_cr_);
    auto coef1 = std::exp(-1.878 * (temperature/T_cr_)) * 0.274/(p_cr_*p_cr_); 
    vector_t Z = 1.0 - coef0 *  pressure.array() 
                        + coef1 * pressure.array().square() ;

    return Z;
}


std::pair<vector_t, vector_t>
papay::speed_of_sound(linearized_fluid_solver *lfs)
{
    auto Z_nodes = compute( lfs->temperature(), lfs->pressure_nodes());
    auto Z_pipes = compute( lfs->temperature(), lfs->pressure_pipes());

    vector_t c2_nodes = Z_nodes.cwiseProduct(R_nodes_) * lfs->temperature(); 
    vector_t c2_pipes = Z_pipes.cwiseProduct(R_pipes_) * lfs->temperature();

    compute_density(lfs,c2_pipes);

    return std::make_pair(c2_nodes, c2_pipes);
}


// ----------------------------------------------------------------------------
//                    GERG equation of state
// ----------------------------------------------------------------------------

gerg_params::gerg_params(const gerg_reducing_params_t&  rp, 
                const gerg_pseudo_critical_pt_t& psc,
                const gerg_thermo_params_t& tp): 
                reducing_params(rp),
                pseudo_critical_pt (psc),
                params (tp) {}; 


gerg_params:: gerg_params(){};


gerg::gerg()
{
    mmi_gerg = vector_t::Zero(21); 

    mmi_gerg <<
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
gerg::compute(  const double& temperature,
                    const vector_t& pressure,
                    const matrix_t& x,
                    const gerg_params& gp)
{
    auto eos =  GERG::thermodynamic_properties(pressure, temperature, x,
                                              gp.reducing_params,
                                              gp.pseudo_critical_pt,
                                              gp.params);
    return eos.Z;
}


std::pair<vector_t, vector_t>
gerg::speed_of_sound(linearized_fluid_solver *lfs)
{
    vector_t Z_nodes = compute( lfs->temperature(),
                                lfs->pressure_nodes(),
                                lfs->x_nodes(), params_nodes_);
    vector_t Z_pipes = compute( lfs->temperature(),
                                lfs->pressure_pipes(),
                                lfs->x_pipes(), params_pipes_);

    vector_t c2_nodes = Z_nodes.cwiseProduct(R_nodes_) * lfs->temperature(); 
    vector_t c2_pipes = Z_pipes.cwiseProduct(R_pipes_) * lfs->temperature();

    compute_density(lfs, c2_pipes);
    return std::make_pair(c2_nodes, c2_pipes);
}


void
gerg::compute_molar_mass(const matrix_t& y_nodes, const matrix_t& y_pipes)
{
    mm_nodes_ = vector_t::Zero(y_nodes.rows()); 
    mm_pipes_ = vector_t::Zero(y_pipes.rows()); 

    std::vector<int> gas_name = 
       {GAS_TYPE::CH4, GAS_TYPE::N2, GAS_TYPE::CO2, GAS_TYPE::C2H6,
        GAS_TYPE::C3H8, GAS_TYPE::i_C4H10, GAS_TYPE::n_C4H10,
        GAS_TYPE::i_C5H12, GAS_TYPE::n_C5H12,GAS_TYPE::C6H14,
        GAS_TYPE::C7H16, GAS_TYPE::C8H18, GAS_TYPE::C9H20,
        GAS_TYPE::C10H22,GAS_TYPE::H2, GAS_TYPE::O2, GAS_TYPE::CO,
        GAS_TYPE::H2O,GAS_TYPE::H2S,GAS_TYPE::He,GAS_TYPE::Ar}; 

    for(size_t i = 0; i <= 20; i++)
        mm_nodes_ +=  mmi_gerg(i) * y_nodes.col(gas_name[i]); 

    for(size_t i = 0; i <= 20; i++)
        mm_pipes_ +=  mmi_gerg(i) * y_pipes.col(gas_name[i]); 

    R_nodes_= Runiversal_ * mm_nodes_.cwiseInverse();     
    R_pipes_= Runiversal_ * mm_pipes_.cwiseInverse();     
}

} //end namespace shimmer