#include <boost/graph/graphviz.hpp>

#include "infrastructure_graph.h"

void write_graphviz(const std::string& fn, const infrastructure_graph& igraph)
{
    std::ofstream ofs("test_graph.dot");
    boost::write_graphviz(ofs, igraph, vertex_property_writer{}, edge_property_writer{});
}