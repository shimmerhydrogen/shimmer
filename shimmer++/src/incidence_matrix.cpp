/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */


#include "incidence_matrix.h"

namespace shimmer{
    
    void
    incidence::compute_triplets(const infrastructure_graph& g)
    {
        auto edge_range = edges(g);
        for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
            auto pipe = g[*itor];   
            auto u = source(*itor, g);
            triplets_in_.push_back(triplet_t(g[u].node_num, pipe.branch_num, 1.0));
 
            auto v = target(*itor, g);
            triplets_out_.push_back(triplet_t(g[v].node_num, pipe.branch_num, 1.0));            
        }

        triplets_ = triplets_out_;
        std::transform(triplets_.cbegin(), triplets_.cend(),
                    triplets_.begin(), [](triplet_t t){    
                        return triplet_t(t.row(), t.col(), -t.value()  );
                    });
       
        auto end = triplets_.end();
        triplets_.insert(end, triplets_in_.begin(),triplets_in_.end());
    }


    void
    incidence::compute_matrix(const infrastructure_graph& g)
    {       
        mat_in_.resize(num_vertices(g), num_edges(g));
        mat_in_.setFromTriplets(triplets_in_.begin(), triplets_in_.end());

        mat_out_.resize(num_vertices(g), num_edges(g));
        mat_out_.setFromTriplets(triplets_out_.begin(), triplets_out_.end()); 

        mat_ = mat_in_ - mat_out_; 
    }


    const sparse_matrix_t& incidence::matrix()      { return mat_;}
    const sparse_matrix_t& incidence::matrix_in()   { return mat_in_;}
    const sparse_matrix_t& incidence::matrix_out()  { return mat_out_;}   
    const sparse_matrix_t& incidence::matrix()      const { return mat_;}
    const sparse_matrix_t& incidence::matrix_in()   const { return mat_in_;}
    const sparse_matrix_t& incidence::matrix_out()  const { return mat_out_;}   


    


} //end namespace shimmer