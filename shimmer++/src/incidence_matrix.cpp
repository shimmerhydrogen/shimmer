/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */


#include "incidence_matrix.h"

namespace shimmer{
    
sparse_matrix_t
incidence_matrix_out(const infrastructure_graph& g)
{
    using triplet_t = Eigen::Triplet<double>;

    std::vector<triplet_t> triplets;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
        auto pipe = g[*itor];   
        auto v = target(*itor, g);
        triplets.push_back(triplet_t(g[v].node_num, pipe.branch_num, double(1)));
    }

    sparse_matrix_t mat(num_vertices(g), num_edges(g));
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    return mat;
}


sparse_matrix_t
incidence_matrix_in(const infrastructure_graph& g)
{
    using triplet_t = Eigen::Triplet<double>;

    std::vector<triplet_t> triplets;
    auto edge_range = edges(g);
    for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
        auto pipe = g[*itor];   
        auto u = source(*itor, g);
        triplets.push_back(triplet_t(g[u].node_num, pipe.branch_num, double(1)));
    }

    sparse_matrix_t mat(num_vertices(g), num_edges(g));
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    return mat;
}


sparse_matrix_t
incidence_matrix(const infrastructure_graph& g)
{
    return  incidence_matrix_in(g) - incidence_matrix_out(g);
}

} //end namespace shimmer