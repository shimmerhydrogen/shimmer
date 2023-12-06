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
    
    pm  = (2.0/3.0) * (io_p - i_p.cwiseProduct(o_p).cwiseQuotient(io_p));
    
    //for(size_t i = 0; i < pm.size(); i++)
    //    std::cout << pm[i] << std::endl;
}


auto
phi_triplets(const double & dt, const double& c2, const infrastructure_graph& g)
{
    std::vector<triplet_t> triplets; 

    int i = 0;
    auto v_range = vertices(g);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        triplets.push_back(triplet_t(i, i, volume(*itor, g) / (c2 * dt)));

    return triplets;
}


void
phi_matrix(const double & dt, const double& c2, const infrastructure_graph& g, sparse_matrix_t& mat)
{
    vector_t phi (num_vertices(g));

    int i = 0;
    auto triplets = phi_triplets(dt, c2, g);   
    mat.resize(num_vertices(g), num_vertices(g));   
    mat.setFromTriplets(triplets.begin(), triplets.end());    
    return;
}

vector_t
compute_expS(const double & c2, const infrastructure_graph& g)
{   
    vector_t  exp_s (num_edges(g)); 
    double gravity = 9.8;
    double factor = 2.0 *  gravity / c2;

    size_t i = 0;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++,i++ )
    {
        auto pipe = g[*itor];   
        auto node_in  = source(*itor, g);
        auto node_out = target(*itor, g);
        auto s = g[node_out].height - g[node_in].height;
        exp_s[i] =  std::exp(factor * s * 0.5 ); 
    }

    return exp_s;       
}


auto 
adp_triplets(const double & c2, const infrastructure_graph& g,
             const incidence& inc)
{
    vector_t exp_s  = compute_expS(c2, g);

    // Referred to pressure p (not squared)

    std::cout << "inc triplets size " << inc.triplets_size() << std::endl;
    std::cout << " exp_s size" << exp_s.size() << std::endl;
     
    auto triplets_out = inc.triplets_out(); 
    auto begin_out = triplets_out.cbegin();
    auto end_out   = triplets_out.cend();

    std::vector<triplet_t> triplets (triplets_out.size());

    std::transform(begin_out, end_out, triplets.begin(),
        [&](triplet_t t){ 
            return triplet_t(t.row(), t.col(), -t.value() * exp_s(t.col())); 
        } 
    );

    auto begin_in = inc.triplets_in().begin();
    auto end_in   = inc.triplets_in().end();
    triplets.insert(triplets.end(), begin_in, end_in);


    // transpose 
    std::transform(triplets.cbegin(), triplets.cend(), triplets.begin(),
        [&](triplet_t t){ 
            return triplet_t(t.col(), t.row(), t.value()); 
        } 
    );


    return triplets;
}


void
adp_matrix(const double & c2, const infrastructure_graph& g,
            const incidence& inc,
            sparse_matrix_t& mat)
{
    // Referred to pressure p (not squared)

    // computation based on triplets
    mat.resize(num_edges(g), num_vertices(g));
    auto triplets = adp_triplets(c2, g, inc);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return;
}

/*
Deprecated since it allows adp_triplets and adp_matrix to be different.
void
adp_matrix( const double & c2, const infrastructure_graph& g,
            const incidence& inc,
            sparse_matrix_t& mat)
{ 
    // Referred to pressure p (not squared)
    auto exp_s  = compute_expS(c2, g);
    
    mat.resize(num_vertices(g), num_edges(g));
    // ==============================================
    // computation based on matrices      
    sparse_matrix_t sE(num_edges(g), num_edges(g));
    sE.setIdentity();
    sE.diagonal() = exp_s;    
    
    mat = (incidence_in - incidence_out * sE).transpose();

    return;
}
*/


auto
resistance_triplets(const double & dt, const double& c2,
                  const vector_t & flux,
                  const vector_t & mean_pressure,
                  const infrastructure_graph  & g)
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

    return triplets;
}


void
resistance_matrix(const double & dt, const double& c2,
                  const vector_t & flux,
                  const vector_t & mean_pressure,
                  const infrastructure_graph  & g,
                  sparse_matrix_t& mat )
{
    auto triplets =  resistance_triplets(dt, c2, flux, mean_pressure, g);
    mat.resize(num_edges(g), num_edges(g));
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    return;
}

} //end namespace shimmer


 // vector_t fromP2p = adp.cwiseAbs().transpose() * pressure;  