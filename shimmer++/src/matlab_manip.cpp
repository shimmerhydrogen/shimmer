 /* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it  
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