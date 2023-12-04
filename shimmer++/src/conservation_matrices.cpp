/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#pragma once

#include <Eigen/Sparse>
#include "infrastructure_graph.h"
#include "../src/geometry_properties.cpp"
#include "../src/incidence_matrix.h"


template<typename T> 
using sparse_matrix_t = Eigen::SparseMatrix<T>; 
template<typename T> 
using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>; 


template<typename T> 
void
average(const  vector_t<T>& pressure, const sparse_matrix_t<T>& incidence_in,
        const sparse_matrix_t<T>& incidence_out,  vector_t<T>& pm)
{
    vector_t<T> i_p = incidence_in.transpose()  * pressure;
    vector_t<T> o_p = incidence_out.transpose()  * pressure;
    vector_t<T> io_p  = i_p + o_p; 
    
    std::cout << "-------------------------------"<< std::endl;
    pm  = (2.0/3.0) * (io_p - i_p.cwiseProduct(o_p).cwiseQuotient(io_p));
    for(size_t i = 0; i < pm.size(); i++)
        std::cout << pm[i] << std::endl;

}


template<typename T> 
void
phi_matrix(const double & dt, const double& c2, const undirected_graph& g, sparse_matrix_t<T>& mat)
{
    auto c2_over_dt = c2/dt;

    std::vector<double> phi (num_edges(g));

    int i = 0;
    auto v_range = vertices(g);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        phi.at(i) = volume(*itor, g) * c2_over_dt;

    mat.setIdentity();
    mat.diagonal() = phi;

    //sparse_matrix_t<T> diag (num_vertices(g), num_vertices(g));    
    //diag.setFromTriplets(triplets.begin(), triplets.end());    
    return;
}


template<typename T> 
void adp_matrix(const T & c2, const undirected_graph& g,
                const sparse_matrix_t<T>& incidence_in,
                const sparse_matrix_t<T>& incidence_out,
                sparse_matrix_t<T>& mat)
{
    // Referred to pressure p (not squared)

    double gravity = 9.8;

    Eigen::Matrix<T, Eigen::Dynamic, 1>  exp_s (num_edges(g)); 

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
       
    sparse_matrix_t<T> sE(num_edges(g), num_edges(g));
    sE.setIdentity();
    sE.diagonal() = exp_s;    
    
    mat = incidence_in - incidence_out * sE;

    return;
}


template <typename T>
using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;  

template<typename T> 
void
resistance_matrix(const T & dt, const T& c2,
                  const vector_t<T> & flux,
                  const vector_t<T> & mean_pressure,
                  const undirected_graph  & g,
                  sparse_matrix_t<T>& mat )
{
    using triplet_t = Eigen::Triplet<T>;
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


 // vector_t<T> fromP2p = adp.cwiseAbs().transpose() * pressure;  