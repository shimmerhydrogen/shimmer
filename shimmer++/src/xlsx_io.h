/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

#pragma once

#include "sol/sol.hpp"
#include "infrastructure_graph.h"

int import_infrastructure_from_xlsx(sol::state&, shimmer::infrastructure_graph&);