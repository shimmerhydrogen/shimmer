 /* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */
#include "../src/matrix_manipulations.h"

namespace shimmer{


std::vector<triplet_t>
build_triplets(const vector_t& vec, size_t row_offset,
                                         size_t col_offset)
{
    std::vector<triplet_t> triplets;

    for(size_t i = 0; i < vec.size(); i++)
        triplets.push_back(triplet_t(i+row_offset, i+col_offset, vec(i)));
    return triplets;    
}


sparse_matrix_t
build_matrix(const vector_t& vec)
{
    sparse_matrix_t mat(vec.size(), vec.size());
    mat.setIdentity();
    mat.diagonal() = vec;   
    return mat;
}


std::vector<triplet_t>
build_triplets(const sparse_matrix_t& sMat, size_t row_offset,
                                            size_t col_offset)
{
    std::vector<triplet_t> triplets;
    size_t count = 0;
    for (int k = 0; k < sMat.outerSize(); ++k)
        for (itor_t it(sMat, k); it; ++it, count++)
            triplets.push_back(triplet_t(it.row() + row_offset,
                                         it.col() + col_offset,
                                         it.value()));
    return triplets;
}

}// end namespace shimmer