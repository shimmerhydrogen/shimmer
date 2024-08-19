/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 *
 * The authors (C) 2023
 */

#pragma once

#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <Eigen/Dense>
#include "../src/boundary.h"
#include "../src/nonpipe_over_edges.h"

namespace shimmer
    {

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

enum GAS_TYPE
{   CH4,
    N2,
    CO2,
    C2H6,
    C3H8,
    i_C4H10,
    n_C4H10,
    i_C5H12,
    n_C5H12,
    C6H14,
    C7H16,
    C8H18,
    C9H20,
    C10H22,
    H2,
    O2,
    CO,
    H2O,
    H2S,
    He,
    Ar,
};

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
    vector_t        gas_mixture;
    std::unique_ptr<station> node_station;

    vertex_properties(){}
    vertex_properties(std::string     iname,
    int             inode_num,
    double          ipressure,
    double          imass_flow,
    double          iheight): name(iname), node_num(inode_num),
                            pressure(ipressure),mass_flow(imass_flow),
                            height(iheight)   {}


    friend std::ostream& operator<<(std::ostream& ofs, const vertex_properties& vp)
    {
        ofs << " name : \"" << vp.name << "\", ";
        ofs << " node_num : " << vp.node_num << ",";
        ofs << " pressure : " << vp.pressure << ", ";
        ofs << " mass_flow: " << vp.mass_flow << "\n";
        if (vp.node_station) {
            ofs << " station:   ";
            vp.node_station->print();
        }

        return ofs;
    }

};

enum class edge_type {
    pipe,
    compressor,
    regulator,
    valve,
};

struct edge_properties {
    edge_type   type;
    int         branch_num;
    double      length;
    double      diameter;
    double      friction_factor;
    int         grid_pts;
    edge_station::station pipe_station;

    double area()   const { return M_PI * 0.25 * diameter * diameter;}
    double volume() const { return area() * length;}

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