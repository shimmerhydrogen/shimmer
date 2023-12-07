/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */


#include "../src/conservation_matrices.h"

namespace shimmer{


void
average(const  vector_t& pressure, const sparse_matrix_t& incidence_in,
        const sparse_matrix_t& incidence_out,  vector_t& pm)
{
    vector_t i_p = incidence_in.transpose()  * pressure;
    vector_t o_p = incidence_out.transpose()  * pressure;
    vector_t io_p  = i_p + o_p; 
    
    std::cout << "-------------------------------"<< std::endl;
    pm  = (2.0/3.0) * (io_p - i_p.cwiseProduct(o_p).cwiseQuotient(io_p));
    for(size_t i = 0; i < pm.size(); i++)
        std::cout << pm[i] << std::endl;

}


void
phi_matrix(const double & dt, const double& c2, const infrastructure_graph& g, sparse_matrix_t& mat)
{
    vector_t phi (num_vertices(g));

    int i = 0;
    auto v_range = vertices(g);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        phi(i) = volume(*itor, g) / (c2 * dt);

    mat.setIdentity();
    mat.diagonal() = phi;

    //sparse_matrix_t diag (num_vertices(g), num_vertices(g));    
    //diag.setFromTriplets(triplets.begin(), triplets.end());    
    return;
}



void 
adp_matrix(const double & c2, const infrastructure_graph& g,
                const sparse_matrix_t& incidence_in,
                const sparse_matrix_t& incidence_out,
                sparse_matrix_t& mat)
{
    // Referred to pressure p (not squared)

    double gravity = 9.8;

    vector_t  exp_s (num_edges(g)); 

    size_t i = 0;
    auto edge_range = edges(g);

    auto factor = 2.0 *  gravity / c2;

    for(auto itor = edge_range.first; itor != edge_range.second;itor++,i++ )
    {
        auto pipe = g[*itor];   
        auto node_in  = source(*itor, g);
        auto node_out = target(*itor, g);
        auto s = g[node_out].height - g[node_in].height;
        exp_s[i] =  std::exp(factor * s * 0.5 ); 
    }
       
    sparse_matrix_t sE(num_edges(g), num_edges(g));
    sE.setIdentity();
    sE.diagonal() = exp_s;    
    
    mat = incidence_in - incidence_out * sE;

    return;
}



void
resistance_matrix(const double & dt, const double& c2,
                  const vector_t & flux,
                  const vector_t & mean_pressure,
                  const infrastructure_graph  & g,
                  sparse_matrix_t& mat )
{
    using triplet_t = Eigen::Triplet<double>;
    std::vector<triplet_t> triplets;

    size_t count = 0;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++,count++ ){
        auto pipe = g[*itor];   
        auto pm = mean_pressure[count];
        auto Omega = 2.0 * pipe.friction_resistance(c2) * std::abs(flux(count)) 
                                + pipe.inertia_resistance(dt,pm); 
        triplets.push_back(triplet_t(count, count, Omega));
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    return;
}

} //end namespace shimmer


 // vector_t fromP2p = adp.cwiseAbs().transpose() * pressure;  