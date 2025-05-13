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
        auto& pipe = g[*itor];
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

