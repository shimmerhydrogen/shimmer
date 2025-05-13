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

#include "solver/variable.h"

namespace shimmer{

variable::variable(){};
variable::variable(const vector_t&p,const vector_t&f,const vector_t&l)
    {
        pressure = p;
        flux = f;
        L_rate = l;
    };


vector_t
variable::make_vector() const
{
    size_t num_pipes = flux.size();
    size_t num_nodes = pressure.size();

    vector_t vec(2 * num_nodes + num_pipes);
    vec.head(num_nodes) = pressure;
    vec.segment(num_nodes, num_pipes) = flux;
    vec.tail(num_nodes) = L_rate;

    return vec;
}

}