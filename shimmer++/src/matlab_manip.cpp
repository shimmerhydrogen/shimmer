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

#include "../src/matlab_manip.h"

namespace shimmer{


matrix_t
build_x_nodes(const infrastructure_graph& g)
{
    size_t num_nodes = num_vertices(g);
    Eigen::MatrixXd x(num_nodes, 21);

    auto index = get(boost::vertex_index, g);

    for (auto vp = vertices(g); vp.first != vp.second; ++vp.first) {
        auto idx = index[*vp.first]; 
        x.row(idx) = g[idx].gas_mixture;
    }

/*
    for(size_t i = 0; i < num_vertices(g); i++) {
        x.row(i) = g[i].gas_mixture;
        std::cout << g[i].node_num <<  " ";
    }
    std::cout << std::endl;
*/
    return x;
}

} // namespace shimmer