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

#include <iostream>
#include <Eigen/Sparse>

namespace shimmer{


using sparse_matrix_t = Eigen::SparseMatrix<double>; 
using triplet_t = Eigen::Triplet<double>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using itor_t = sparse_matrix_t::InnerIterator;


std::vector<triplet_t>
build_triplets(const vector_t& vec, size_t row_offset = 0,
                                    size_t col_offset = 0);


sparse_matrix_t
build_matrix(const vector_t& vec);


std::vector<triplet_t>
build_triplets(const sparse_matrix_t& sMat, size_t row_offset = 0,
                                            size_t col_offset = 0);

}// end namespace shimmer
