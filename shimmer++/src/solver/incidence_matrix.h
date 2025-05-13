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

#pragma once

#include <Eigen/Sparse>
#include "infrastructure_graph.h"

namespace shimmer{

using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;


class incidence
{
    sparse_matrix_t mat_; 
    sparse_matrix_t mat_in_;
    sparse_matrix_t mat_out_; 

    std::vector<triplet_t> triplets_;
    std::vector<triplet_t> triplets_in_;
    std::vector<triplet_t> triplets_out_;


    void
    compute_triplets(const infrastructure_graph& g);

    void
    compute_matrix(const infrastructure_graph& g);

public:
    incidence(const infrastructure_graph& g)
    {
        compute_triplets(g);
        compute_matrix(g);
    }

    incidence(){}

    const sparse_matrix_t& matrix();     
    const sparse_matrix_t& matrix_in();   
    const sparse_matrix_t& matrix_out();   
    const sparse_matrix_t& matrix() const;      
    const sparse_matrix_t& matrix_in()  const;   
    const sparse_matrix_t& matrix_out() const;  

};



} //end namespace shimmer