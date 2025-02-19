/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#include <cassert>
#include "../src/conservation_matrices.h"
#include "../src/matrix_manipulations.h"

namespace shimmer{


vector_t
average(const vector_t& pressure, const incidence& inc)
{
    vector_t i_p = inc.matrix_in().transpose()  * pressure;
    vector_t o_p = inc.matrix_out().transpose() * pressure;
    vector_t io_p  = i_p + o_p; 
    
    return (2.0/3.0) * (io_p - i_p.cwiseProduct(o_p).cwiseQuotient(io_p));

}


vector_t
phi_vector(const double& dt, const vector_t& c2, const infrastructure_graph& g)
{
    vector_t phi = vector_t::Zero(num_vertices(g));

    int i = 0;
    auto v_range = vertices(g);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        phi(i) = volume(*itor, g) / (c2(i) * dt);
    
    return phi;
}


sparse_matrix_t
phi_matrix(const double& dt, const vector_t& c2, const infrastructure_graph& g)
{
    vector_t phi_vec = phi_vector(dt, c2, g);   

    return build_matrix(phi_vec);
}


vector_t
compute_expS(const vector_t & c2, const infrastructure_graph& g)
{   
    vector_t  exp_s (num_edges(g)); 
    double gravity = 9.8;

    size_t i = 0;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++,i++)
    {
        auto node_in  = source(*itor, g);
        auto node_out = target(*itor, g);
        auto s = g[node_out].height - g[node_in].height;
        exp_s(i) =  std::exp(s * gravity / c2(i)); 
    }

    return exp_s;       
}


sparse_matrix_t
adp_matrix( const vector_t & c2, const infrastructure_graph& g,
            const incidence& inc)
{ 
    // Referred to pressure p (not squared)
    vector_t exp_s  = compute_expS(c2, g);
    sparse_matrix_t sE = build_matrix(exp_s);  

    sparse_matrix_t ADP(num_edges(g), num_vertices(g));
    ADP = (inc.matrix_in() - inc.matrix_out() * sE).transpose();

    return ADP;
}


vector_t
resistance_inertia( const double & dt, const vector_t & pipes_pressure,
                    const incidence& inc, const infrastructure_graph  & g)
{
    vector_t Omega = vector_t::Zero(num_edges(g));
    assert(pipes_pressure.size() == num_edges(g));

    size_t i = 0;
    auto edge_range = edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;
    for(auto itor = begin; itor != end; itor++,i++ ){
        const auto& pipe = g[*itor];   
        auto pm = pipes_pressure(i);
        Omega(i) = inertia_resistance(pipe, dt, pm); 
    }
    return Omega;
} 


vector_t
resistance_friction(const double& temperature,
                    const vector_t& mu,
                    const vector_t& c2,
                    const vector_t & flux,
                    const infrastructure_graph  & g)
{
    vector_t Omega = vector_t::Zero(num_edges(g));

    auto edge_range = edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;

    size_t i = 0;
    for(auto itor = begin; itor != end; itor++,i++ ){
        const auto& pipe = g[*itor]; 
        auto rf = friction_resistance(pipe, temperature, mu(i),c2(i), flux(i));       
        Omega(i) = 2.0 * rf * std::abs(flux(i)); 
    }
    return Omega;
} 


} //end namespace shimmer


 // vector_t fromP2p = adp.cwiseAbs().transpose() * pressure;  