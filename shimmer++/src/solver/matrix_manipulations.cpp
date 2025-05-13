/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "solver/matrix_manipulations.h"

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