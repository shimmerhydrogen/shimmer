/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2024
 * karol.cascavita@polito.it  
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "../src/infrastructure_graph.h"
#include "../src/boundary.h"
#include "verify_test.h"

using namespace shimmer;

size_t num_steps = 25;
size_t num_inlet = 1;
size_t num_pipes = 15;
size_t num_nodes = 13;

/*
void
make_init_graph(infrastructure_graph& g, const  std::vector<std::unique_ptr<station>>& stations)
{

    std::vector<vertex_descriptor> vds;

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( boost::add_vert = 0( { "station 10",10, 66.569303493,  15 ,0.,x , stations[10]}, g) );
    vds.push_back( boost::add_vertex( { "station 11",11, 66.569303493,  40 ,0.,x , stations[11]}, g) );
    vds.push_back( boost::add_vertex( { "station 12",12, 67.241720700,  10 ,0.,x , stations[12]}, g) );


    using eprop_t = edge_properties;

    edge_properties ep0  = {edge_type::pipe, 0,  80000,	1.2,	1.20E-05};
    edge_properties ep1  = {edge_type::pipe, 1,  16000,	0.6,	1.20E-05};
    edge_properties ep2  = {edge_type::pipe, 2,  40000,	0.8,	1.20E-05};
    edge_properties ep3  = {edge_type::pipe, 3, 160000,	0.7,	1.20E-05};
    edge_properties ep4  = {edge_type::pipe, 4, 200000,	0.8,	1.20E-05};
    edge_properties ep5  = {edge_type::pipe, 5,  24000,	0.6,	1.20E-05};
    edge_properties ep6  = {edge_type::pipe, 6, 120000,	0.2,	1.20E-05};
    edge_properties ep7  = {edge_type::pipe, 7,  80000,	0.9,	1.20E-05};
    edge_properties ep8  = {edge_type::pipe, 8,  64000,	0.7,	1.20E-05};
    edge_properties ep9  = {edge_type::pipe, 9, 240000,	0.6,	1.20E-05};
    edge_properties ep10 = {edge_type::pipe,10,  28000,	0.2,	1.20E-05};
    edge_properties ep11 = {edge_type::pipe,11,  80000,	0.9,	1.20E-05};
    edge_properties ep12 = {edge_type::pipe,12, 160000,	0.7,	1.20E-05};
    edge_properties ep13 = {edge_type::pipe,13,  40000,	0.3,	1.20E-05};
    edge_properties ep14 = {edge_type::pipe,14, 320000,	0.9,	1.20E-05};


    boost::add_edge( vds[ 0], vds[ 3], ep0, g);
    boost::add_edge( vds[ 1], vds[ 2], ep1, g);
    boost::add_edge( vds[ 2], vds[ 3], ep2, g);
    boost::add_edge( vds[ 2], vds[ 4], ep3, g);
    boost::add_edge( vds[ 3], vds[ 4], ep4, g);
    boost::add_edge( vds[ 4], vds[ 5], ep5, g);
    boost::add_edge( vds[ 4], vds[ 7], ep6, g);
    boost::add_edge( vds[ 6], vds[ 4], ep7, g);
    boost::add_edge( vds[ 7], vds[ 6], ep8, g);
    boost::add_edge( vds[11], vds[ 6], ep9, g);
    boost::add_edge( vds[12], vds[ 7], ep10, g);
    boost::add_edge( vds[ 8], vds[ 7], ep11, g);
    boost::add_edge( vds[ 7], vds[ 9], ep12, g);
    boost::add_edge( vds[ 9], vds[10], ep13, g);
    boost::add_edge( vds[ 3], vds[ 9], ep14, g);

}
*/



enum station_type
{
    INLET,
    OUTLET,
    JUNCTION,
};


std::vector<station_type>
make_stations_type_vector()
{
    vector_t outlet_nodes(9);
    outlet_nodes << 1, 2, 5, 6, 8, 9, 10, 11, 12;

    vector_t inlet_nodes(1);
    inlet_nodes << 0; 

    vector_t junction_nodes(3);
    junction_nodes << 3, 4, 7; 

    std::vector<station_type> vec(num_nodes); 
    for(size_t i = 0; i < outlet_nodes.size(); i++)
        vec[outlet_nodes[i]] = station_type::OUTLET;            
    for(size_t i = 0; i < inlet_nodes.size(); i++)
        vec[inlet_nodes[i]] = station_type::INLET;            
    for(size_t i = 0; i < junction_nodes.size(); i++)
        vec[junction_nodes[i]] = station_type::JUNCTION;            
    
    return vec;
}


std::pair<vector_t, matrix_t>
read_boundary_data()
{
    vector_t Pset(num_steps);
    Pset << 70.000000000000000,  69.416666666666671,  68.833333333333343,
            68.250000000000000,  67.666666666666671,  67.083333333333343,
            66.500000000000000,  67.083333333333329,  67.666666666666657,
            68.249999999999986,  68.833333333333329,  69.416666666666657,
            70.000000000000000,  72.333333333333329,  74.666666666666657,
            76.999999999999972,  79.333333333333300,  81.666666666666629,
            84.000000000000000,  81.666666666666657,  79.333333333333329,
            77.000000000000000,  74.666666666666686,  72.333333333333357,
            70.000000000000000;
    Pset *=1E5;
    //---------------------------------------------------------------
    matrix_t Gsnam(num_steps, num_nodes);
    std::ifstream ifs("../unit_tests/gsnam.txt");
    if (!ifs.is_open()) 
    {
        std::cout<< "ERROR: gsnam.txt not open" << std::endl;
        exit(1);
    }

    for (size_t icol = 0; icol < num_nodes; icol++)
        for (size_t irow = 0; irow < num_steps; irow++)
            ifs >>  Gsnam(irow, icol);
    ifs.close();

    return std::make_pair(Pset, Gsnam);
}
int main()
{
    std::vector<double> ref_values = {7000000.0,  25.0, 25.0, 0.0,
                                      0.0,  20.0, 30.0,  0.0, 62.5,
                                      20.0, 16.5, 30.0, 10.0};
    //---------------------------------------------------------------
    auto [Pset, Gsnam] = read_boundary_data();
    //---------------------------------------------------------------
    auto station_type_vec = make_stations_type_vector();
    //---------------------------------------------------------------
    std::vector<std::unique_ptr<station>> stations(num_nodes);

    for(size_t i = 0 ; i < num_nodes; i++)
    {   
        switch(station_type_vec[i])
        {
            case(station_type::INLET):{
                auto s = make_inlet(Pset);
                stations[i] = std::make_unique<inlet_station>();
                stations[i]->set_state(s);
                break;
            }
            case(station_type::OUTLET):{
                auto s = make_outlet(Gsnam.col(i));
                stations[i] = std::make_unique<outlet_station>();
                stations[i]->set_state(s);
                break;
            }
            case(station_type::JUNCTION):
                stations[i] = std::make_unique<junction>();
                break;
            default:
                throw std::invalid_argument("Station type not found");    
        }
    }
    //---------------------------------------------------------------

    infrastructure_graph g;
    std::vector<vertex_descriptor> vds;

    auto add_vertex = [&](vertex_properties&& vp, const vector_t& x_in, size_t i) 
    {
        vp.gas_mixture = x_in;
        vp.node_station = std::move(stations[i]);
        auto v = boost::add_vertex(g);
        g[v] = std::move(vp);
        return v;
    };

    vector_t  x = vector_t::Zero(21);
    x(GAS_TYPE::CH4) = 1.0;

    vds.push_back( add_vertex(vertex_properties("station 0",  0, 70.000000000,-230 ,0.),x ,0));
    vds.push_back( add_vertex(vertex_properties("station 1",  1, 70.000000000,  20 ,0.),x ,1));
    vds.push_back( add_vertex(vertex_properties("station 2",  2, 69.300000000,  25 ,0.),x ,2));
    vds.push_back( add_vertex(vertex_properties("station 3",  3, 69.300000000,   0 ,0.),x ,3));
    vds.push_back( add_vertex(vertex_properties("station 4",  4, 68.607000000,   0 ,0.),x ,4));
    vds.push_back( add_vertex(vertex_properties("station 5",  5, 67.920930000,  20 ,0.),x ,5));
    vds.push_back( add_vertex(vertex_properties("station 6",  6, 67.241720700,  30 ,0.),x ,6));
    vds.push_back( add_vertex(vertex_properties("station 7",  7, 67.920930000,   0 ,0.),x ,7));
    vds.push_back( add_vertex(vertex_properties("station 8",  8, 67.241720700,  50 ,0.),x ,8));
    vds.push_back( add_vertex(vertex_properties("station 9",  9, 67.241720700,  20 ,0.),x ,9));
    vds.push_back( add_vertex(vertex_properties("station 10",10, 66.569303493,  15 ,0.),x ,10));
    vds.push_back( add_vertex(vertex_properties("station 11",11, 66.569303493,  40 ,0.),x ,11));
    vds.push_back( add_vertex(vertex_properties("station 12",12, 67.241720700,  10 ,0.),x ,12));

    using eprop_t = edge_properties;

    edge_properties ep0  = {edge_type::pipe, 0,  80000,	1.2,	1.20E-05};
    edge_properties ep1  = {edge_type::pipe, 1,  16000,	0.6,	1.20E-05};
    edge_properties ep2  = {edge_type::pipe, 2,  40000,	0.8,	1.20E-05};
    edge_properties ep3  = {edge_type::pipe, 3, 160000,	0.7,	1.20E-05};
    edge_properties ep4  = {edge_type::pipe, 4, 200000,	0.8,	1.20E-05};
    edge_properties ep5  = {edge_type::pipe, 5,  24000,	0.6,	1.20E-05};
    edge_properties ep6  = {edge_type::pipe, 6, 120000,	0.2,	1.20E-05};
    edge_properties ep7  = {edge_type::pipe, 7,  80000,	0.9,	1.20E-05};
    edge_properties ep8  = {edge_type::pipe, 8,  64000,	0.7,	1.20E-05};
    edge_properties ep9  = {edge_type::pipe, 9, 240000,	0.6,	1.20E-05};
    edge_properties ep10 = {edge_type::pipe,10,  28000,	0.2,	1.20E-05};
    edge_properties ep11 = {edge_type::pipe,11,  80000,	0.9,	1.20E-05};
    edge_properties ep12 = {edge_type::pipe,12, 160000,	0.7,	1.20E-05};
    edge_properties ep13 = {edge_type::pipe,13,  40000,	0.3,	1.20E-05};
    edge_properties ep14 = {edge_type::pipe,14, 320000,	0.9,	1.20E-05};


    boost::add_edge( vds[ 0], vds[ 3], ep0, g);
    boost::add_edge( vds[ 1], vds[ 2], ep1, g);
    boost::add_edge( vds[ 2], vds[ 3], ep2, g);
    boost::add_edge( vds[ 2], vds[ 4], ep3, g);
    boost::add_edge( vds[ 3], vds[ 4], ep4, g);
    boost::add_edge( vds[ 4], vds[ 5], ep5, g);
    boost::add_edge( vds[ 4], vds[ 7], ep6, g);
    boost::add_edge( vds[ 6], vds[ 4], ep7, g);
    boost::add_edge( vds[ 7], vds[ 6], ep8, g);
    boost::add_edge( vds[11], vds[ 6], ep9, g);
    boost::add_edge( vds[12], vds[ 7], ep10, g);
    boost::add_edge( vds[ 8], vds[ 7], ep11, g);
    boost::add_edge( vds[ 7], vds[ 9], ep12, g);
    boost::add_edge( vds[ 9], vds[10], ep13, g);
    boost::add_edge( vds[ 3], vds[ 9], ep14, g);


    auto v_range = vertices(g);
    for(auto itor = v_range.first; itor != v_range.second; itor++)
        std::cout<<g[*itor] << std::endl;

    size_t i = 0;
    size_t step = 0;

    vector_t values = vector_t::Zero(num_nodes);
    for(auto itor = v_range.first; itor != v_range.second; itor++, i++)
        values[i] =  g[*itor].node_station->boundary().value(step);
    //make_init_graph(graph, node_stations);
    bool pass = verify_test("stations", values, ref_values); 

    std::cout<< !pass << std::endl;
    return !pass;
}
