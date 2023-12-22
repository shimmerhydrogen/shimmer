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
        auto pipe = g[*itor];   
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


sparse_matrix_t
apa_matrix( const vector_t & c2, const  vector_t& pressure, 
            const infrastructure_graph& g, const incidence& inc)
{
    sparse_matrix_t ADP = adp_matrix(c2, g, inc);
    sparse_matrix_t diag_ADP_p = build_matrix( ADP.cwiseAbs() * pressure);
    
    sparse_matrix_t sAPA(num_edges(g), num_vertices(g));
    sAPA = diag_ADP_p * ADP; 

    return sAPA;
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
        auto pipe = g[*itor];   
        auto pm = pipes_pressure(i);
        Omega(i) = inertia_resistance(pipe, dt, pm); 
    }
    return Omega;
} 


vector_t
resistance_friction(const double& temperature, const vector_t& c2,
                    const vector_t & flux,
                    const infrastructure_graph  & g)
{
    vector_t Omega = vector_t::Zero(num_edges(g));

    size_t i = 0;
    auto edge_range = edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;
    for(auto itor = begin; itor != end; itor++,i++ ){
        auto pipe = g[*itor]; 
        auto node_in = g[source(*itor, g)];                        
        auto rf = friction_resistance(pipe, node_in, c2(i), temperature, flux(i));
        Omega(i) = 2.0 * rf * std::abs(flux(i)); 

    }
    return Omega;
} 

/* 
auto
momemtum(const double& dt, const double& temperature,
         const vector_t& flux, const vector_t& flux_old,
         const vector_t& pressure, const incidence& inc,
         const infrastructure_graph & graph)
{
    
    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);
    //size_t num_pipes_ext = num_pipes;

    vector_t pipes_pressure = average(pressure);
    vector_t c2 = speed_of_sound(temperature,  pm);        

    sparse_matrix_t sAPA  = apa_matrix(c2, pressure, graph, inc);

    vector_t rf = resistance_friction(temperature, c2, flux, graph);
    vector_t ri = resistance_inertia(dt, pipes_pressure, inc, graph);

    auto t_sR   = build_triplets(-rf-ri, num_nodes, num_nodes);
    auto t_sAPA = build_triplets( sAPA,  num_nodes, 0);

    std::vector<triplet_t> triplets =  t_sAPA; 
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());
    //triplets.insert(triplets.begin(), t_sIc.begin(), t_sIc.end());

    vector_t rhs_momentum =  -0.5 * rf.array() * flux.array()
                                   - ri.array() * flux_old.array();

    return std::make_pair(triplets, rhs);
}


auto
continuity(const double& dt, const double& temperature,
        const  vector_t& pressure, const vector_t& pressure_old,
        const incidence& inc,
        const infrastructure_graph & graph)
{
    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);

    vector_t c2 = speed_of_sound(temperature,  pressure);   

    vector_t phi_vec = phi_vector(dt, c2, graph);
    auto t_sPHI = build_triplets( phi_vec);
    auto t_sA   = build_triplets( inc.matrix(), 0, num_nodes );

    std::vector<triplet_t> triplets =  t_sPHI; 
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());

    vector_t rhs = phi_vec.array() * pressure_old.array();

    return std::make_pair(triplets, rhs);
}
*/

} //end namespace shimmer


 // vector_t fromP2p = adp.cwiseAbs().transpose() * pressure;  