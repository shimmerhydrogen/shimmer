/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#include <boost/graph/graphviz.hpp>

#include "infrastructure_graph.h"

namespace shimmer{ 
    
void write_graphviz(const std::string& fn, const infrastructure_graph& igraph)
{
    std::ofstream ofs("test_graph.dot");
    boost::write_graphviz(ofs, igraph, vertex_property_writer{}, edge_property_writer{});
}

} //end namespace shimmer