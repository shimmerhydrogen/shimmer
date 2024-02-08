/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */


#include "geometry_properties.h"

namespace shimmer{


double area(const edge_properties& ep)
{
    return ep.area();
}

double volume(const edge_properties& ep)
{
    return ep.volume();
}

vector_t area( const infrastructure_graph& g)
{
    vector_t a(num_edges(g));
    auto edge_range = boost::edges(g);
    auto begin = edge_range.first;
    auto end = edge_range.second;
    size_t i = 0;
    for(auto itor = begin; itor != end; itor++,i++ ){
        auto pipe = g[*itor];
        a(i) = pipe.area();   
    }
    return a;
}

double volume(const infrastructure_graph::vertex_descriptor&  v, const infrastructure_graph& g)
{

    double vol = 0;
    auto edge_out = out_edges(v, g);

    for(auto itor = edge_out.first; itor != edge_out.second;itor++ )
        vol += volume(g[*itor]);

    return vol * 0.5;   
}

} //end namespace shimmer

