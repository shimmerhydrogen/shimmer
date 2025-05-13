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

#include "verify_test.h"
#include <iomanip>

bool verify_test(const std::string & name, 
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& vals,
                 const std::vector<double>& ref )
{
    if(vals.size() != ref.size()) 
        throw std::invalid_argument("ERROR: solution and reference dont have same size.");


    bool pass = true;

   
    for (int k = 0; k < vals.size(); ++k)
    {
        std::cout << std::setprecision(16) << vals[k]  <<  std::endl ;
        double e_val;
        if(std::abs(ref.at(k)) > 1.e-12)
            e_val = std::abs((vals[k] - ref.at(k))) /  std::abs(ref.at(k));
        else 
            e_val = vals[k];
        if(e_val > 1.e-3)
        {
            pass = false;
            break;  
        }
    }
    
    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << "  Test " << name << ".........." <<  passfail(pass) << std::endl;

    if(!pass) 
    {
        std::cout << " * Vec:" << std::endl;
        std::cout << vals << std::endl;
    }

    return pass;
}



bool verify_test(const std::string & name, 
                 const Eigen::SparseMatrix<double>& mat,
                 const std::vector<std::array<double, 3>>& ref )
{

    std::cout << " Inside verify test 2 " << std::endl; 

    if(mat.nonZeros() != ref.size()) 
        throw std::invalid_argument("ERROR: solution and reference dont have same size.");


    using itor_t = Eigen::SparseMatrix<double>::InnerIterator;
    bool pass = true;
    size_t count = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (itor_t it(mat,k); it; ++it, count++)
        { 
            std::cout << std::setprecision(16) << "(" << it.row() 
                      << " , " << it.col() << " , " << it.value() 
                      << " ) " << std::endl ;

            auto t = ref.at(count);
            auto e_val = std::abs((it.value() - t[2])/t[2]);
            if((it.row() != t[0])  || (it.col() != t[1]) || (e_val > 1.e-12))
            {
                pass = false;
                break;  
            }
        }
    }
    
    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << "  Test " << name << " .........." <<  passfail(pass) << std::endl;

    return pass;
}

