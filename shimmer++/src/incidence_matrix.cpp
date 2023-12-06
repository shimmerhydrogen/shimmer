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
    incidence::incidence_in_triplets(const infrastructure_graph& g)
    {
        auto edge_range = edges(g);
        for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
            auto pipe = g[*itor];   
            auto u = source(*itor, g);
            triplets_in_.push_back(triplet_t(g[u].node_num, pipe.branch_num, 1.0));
        }
    }


    void
    incidence::incidence_out_triplets(const infrastructure_graph& g)
    {
        auto edge_range = edges(g);
        for(auto itor = edge_range.first; itor != edge_range.second;itor++ ){
            auto pipe = g[*itor];   
            auto v = target(*itor, g);
            triplets_out_.push_back(triplet_t(g[v].node_num, pipe.branch_num, 1.0));
        }
    }


    void
    incidence::incidence_triplets(const infrastructure_graph& g)
    {
        incidence_in_triplets(g);
        incidence_out_triplets(g);

        triplets_ = triplets_out_;
        std::transform(triplets_.cbegin(), triplets_.cend(),
                    triplets_.begin(), [](triplet_t t){    
                        return triplet_t(t.row(), t.col(), -t.value()  );
                    });
       
        auto end = triplets_.end();
        triplets_.insert(end, triplets_in_.begin(),triplets_in_.end());
    }


    void 
    incidence::incidence_matrix_in(const infrastructure_graph& g)
    {
        mat_in_.resize(num_vertices(g), num_edges(g));
        mat_in_.setFromTriplets(triplets_in_.begin(), triplets_in_.end());
    }


    void
    incidence::incidence_matrix_out(const infrastructure_graph& g)
    {
        mat_out_.resize(num_vertices(g), num_edges(g));
        mat_out_.setFromTriplets(triplets_out_.begin(), triplets_out_.end()); 
    }


    void
    incidence::incidence_matrix(const infrastructure_graph& g)
    {       
        incidence_matrix_in(g);
        incidence_matrix_out(g);

        mat_ = mat_in_ - mat_out_; 
    }



    size_t incidence::triplets_size()    {return triplets_.size();}
    size_t incidence::triplets_in_size() {return triplets_in_.size();}
    size_t incidence::triplets_out_size(){return triplets_out_.size();}
    size_t incidence::triplets_size()    const {return triplets_.size();}
    size_t incidence::triplets_in_size() const {return triplets_in_.size();}
    size_t incidence::triplets_out_size()const {return triplets_out_.size();}

    const sparse_matrix_t& incidence::matrix()      { return mat_;}
    const sparse_matrix_t& incidence::matrix_in()   { return mat_in_;}
    const sparse_matrix_t& incidence::matrix_out()  { return mat_out_;}   
    const sparse_matrix_t& incidence::matrix()      const { return mat_;}
    const sparse_matrix_t& incidence::matrix_in()   const { return mat_in_;}
    const sparse_matrix_t& incidence::matrix_out()  const { return mat_out_;}   

    const std::vector<triplet_t>& incidence::triplets()     const {return triplets_;}
    const std::vector<triplet_t>& incidence::triplets_in()  const  {return triplets_in_;}
    const std::vector<triplet_t>& incidence::triplets_out() const  {return triplets_out_;}

    const std::vector<triplet_t>& incidence::triplets()      {return triplets_;}
    const std::vector<triplet_t>& incidence::triplets_in()    {return triplets_in_;}
    const std::vector<triplet_t>& incidence::triplets_out()   {return triplets_out_;}

    


} //end namespace shimmer