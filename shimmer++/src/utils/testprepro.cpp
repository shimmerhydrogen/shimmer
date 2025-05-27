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

#include "infra/infrastructure.h"
#include "errors.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Please specify database file name" << std::endl;
        return 1;
    }

    shimmer::infrastructure infra;

    int err = shimmer::load(argv[1], infra);
    if (err != SHIMMER_SUCCESS) {
        std::cout << "Problem detected while loading DB" << std::endl;
        return 1;
    }

    shimmer::infrastructure infraqt;
    preprocess_for_quality_tracking(infra, infraqt, 5000);

    return 0;
}