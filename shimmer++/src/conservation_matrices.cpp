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
resistance_friction(const double& temperature,
                    const vector_t& c2,
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


pair_trip_vec_t
continuity( const double& dt, 
            const double& temperature,
            const vector_t& pressure,
            const vector_t& pressure_old,
            const incidence& inc,
            const infrastructure_graph & graph,
            const matrix_t& x,
            const vector_t& RR,
            const gerg_params& gerg)
{
    size_t num_nodes = num_vertices(graph); 
    size_t num_pipes = num_edges(graph);

    auto eos = equation_of_state(temperature, pressure, x, gerg);
    vector_t c2 = eos.Z.cwiseProduct(RR) * temperature; 
 
    vector_t phi_vec = phi_vector(dt, c2, graph);
    auto t_sPHI = build_triplets( phi_vec);
    auto t_sA   = build_triplets( inc.matrix(), 0, num_nodes );

    std::vector<triplet_t> triplets =  t_sPHI; 
    triplets.insert(triplets.begin(), t_sA.begin(), t_sA.end());

    vector_t rhs = phi_vec.array() * pressure_old.array();

    return std::make_pair(triplets, rhs);
}


std::pair<pair_trip_vec_t, vector_t>
momentum(const double& dt,
         const double& temperature,
         const vector_t& flux, 
         const vector_t& flux_old,
         const vector_t& nodes_pressure,
         const incidence& inc,
         const infrastructure_graph & graph,
         const matrix_t& x,
         const vector_t& RR,
         const vector_t& molar_mass,
         const gerg_params& gerg)
{  
    double factor = 1.;//1.e-9;
    size_t num_nodes = num_vertices(graph);
    size_t num_pipes = num_edges(graph);
    //size_t num_pipes_ext = num_pipes;

    vector_t pipes_pressure = average(nodes_pressure, inc);
    auto eos = equation_of_state(temperature, pipes_pressure, x, gerg);
    vector_t c2 = eos.Z.cwiseProduct(RR) * temperature;

    sparse_matrix_t sADP = factor * adp_matrix(c2, graph, inc);
    vector_t ADP_p = sADP.cwiseAbs() * nodes_pressure;

    vector_t rf = resistance_friction(temperature, c2, flux, graph);
    vector_t ri = resistance_inertia(dt, pipes_pressure, inc, graph);
    vector_t r = (-rf-ri).array() / ADP_p.array();
    auto t_sR   = build_triplets( r, num_nodes, num_nodes);
    auto t_sADP = build_triplets( sADP,  num_nodes, 0);

    //std::vector<triplet_t> triplets =  t_sAPA;
    std::vector<triplet_t> triplets =  t_sADP;
    triplets.insert(triplets.begin(), t_sR.begin(), t_sR.end());

    vector_t rhs = factor * ((-0.5) * rf.array() * flux.array()
                                   - ri.array() * flux_old.array())/ADP_p.array();

    vector_t area(num_pipes);
    auto edge_range = boost::edges(graph);
    auto begin = edge_range.first;
    auto end = edge_range.second;
    size_t i = 0;
    for(auto itor = begin; itor != end; itor++,i++ ){
        auto pipe = graph[*itor];
        area(i) = pipe.area();   
    }

    /// rho [kg/m3] actual density of the gas (pipeline based)
    /// vel [m/s] velocity of the gas within pipes.
    vector_t rho = eos.D.cwiseProduct(inc.matrix_in().transpose() * molar_mass); 
    vector_t vel = flux.cwiseQuotient(area.cwiseProduct(rho));

    return std::make_pair(std::make_pair(triplets, rhs), vel);
}


pair_trip_vec_t
boundary(const double& p_in,
        const vector_t& vel,
        const vector_t& flux_ext,
        const incidence& inc,
        const infrastructure_graph& graph)
{
    size_t num_pipes = num_edges(graph);
    size_t num_nodes = num_vertices(graph);

    sparse_matrix_t sId (num_nodes, num_nodes);
    sId.setIdentity();
    auto triplets = build_triplets(sId, 0, num_nodes + num_pipes);

    vector_t rhs = flux_ext;

    if (vel(0) > 0)
    {
        triplets.push_back(triplet_t(num_pipes + num_nodes, 0, 1.));
        sId.coeffRef(0, 0) = 0.0;
        rhs(0) = p_in;
    }
    else
    {
        sId.coeffRef(0, 0) = 0.0; // Warning: This is not the same MATRIX_k(dimn+dimb+1,1)=0; 
        sId.coeffRef(0, 0) = 1.0;
        rhs(0) = 0.0;
    }

    auto t_sId_bcnd = build_triplets(sId, num_nodes + num_pipes, num_nodes + num_pipes);
    triplets.insert(triplets.begin(), t_sId_bcnd.begin(), t_sId_bcnd.end());

    return  std::make_pair(triplets, rhs); 
}


} //end namespace shimmer


 // vector_t fromP2p = adp.cwiseAbs().transpose() * pressure;  