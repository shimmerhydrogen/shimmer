/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#pragma once

#include <boost/graph/adjacency_list.hpp>

namespace shimmer {

struct gas_descriptor {
    std::string     name;
    double          percentage;
};

struct vertex_properties {
    std::string     name;
    int             node_num;
    double          pressure;
    double          mass_flow;
    double          height;
    std::vector<gas_descriptor> gas_mixture;

   
};

enum class edge_type {
    pipe,
    resistor,
    compressor,
    regulator,
    valve
};

struct edge_properties {
    edge_type   type;
    int         branch_num;
    double      length;
    double      diameter;
    double      friction_factor;
    int         grid_pts;

    double area()   const { return M_PI * 0.25 * diameter * diameter;}
    double volume() const { return area() * length;}

    double inertia_resistance(const double& dt, const double& mean_pressure) const
    {
        return  2.0 * length * mean_pressure / (dt * area()); 
    }

    double friction_resistance(const double & c2) const
    {
        double a = area();

        return  friction_factor * c2 * length / (a * a * diameter );    
    }

    friend std::ostream& operator<<(std::ostream& ofs, const edge_properties& ep) {
        ofs << " branch_num : " << ep.branch_num << "\n";
        ofs << " length     : " << ep.length << "\n";
        ofs << " diameter   : " << ep.diameter << "\n";
        ofs << " friction factor : " << ep.friction_factor << "\n";

        return ofs;
    }
};

using infrastructure_graph = boost::adjacency_list<boost::listS,
    boost::vecS, boost::undirectedS, vertex_properties, edge_properties>;

using vertex_descriptor = infrastructure_graph::vertex_descriptor;
using edge_descriptor   = infrastructure_graph::edge_descriptor;

struct vertex_property_writer {
    template<typename Vertex>
    void operator()(std::ostream& ofs, const Vertex& v) const {
        ofs << "[label = \"" << v << "\"]";
    }
};

struct edge_property_writer {
    template<typename Edge>
    void operator()(std::ostream& ofs, const Edge& e) const {
        ofs << "[label = \"" << e << "\"]";
    }
};

void write_graphviz(const std::string&, const infrastructure_graph&);

} //end namespace shimmer