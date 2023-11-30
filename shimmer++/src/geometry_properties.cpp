/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

 #pragma once

#include "infrastructure_graph.h"



double area(const edge_properties& ep)
{
    return M_PI * ep.diameter * ep.diameter * 0.25;
}

double volume(const edge_properties& ep)
{
    return area(ep) * ep.length;
}

double volume(const vertex_descriptor&  v, const infrastructure_graph& g)
{

    double vol = 0;
    auto edge_out = out_edges(v, g);
    //auto edge_in  = in_edges(v, g);

    for(auto itor = edge_out.first; itor != edge_out.second;itor++ )
        vol += volume(g[*itor]);

    //for(auto itor = edge_in.first; itor != edge_in.second;itor++ )
    //    vol += volume(g[*itor]);

    return vol * 0.5;   
}


double volume(const undirected_graph::vertex_descriptor&  v, const undirected_graph& g)
{

    double vol = 0;
    auto edge_out = out_edges(v, g);

    for(auto itor = edge_out.first; itor != edge_out.second;itor++ )
        vol += volume(g[*itor]);

    return vol * 0.5;   
}

