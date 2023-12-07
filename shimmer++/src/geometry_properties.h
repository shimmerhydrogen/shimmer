/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */
# pragma once

#include "infrastructure_graph.h"

double area(const edge_properties& ep);

double volume(const edge_properties& ep);

double volume(const infrastructure_graph::vertex_descriptor&  v, const infrastructure_graph& g);
