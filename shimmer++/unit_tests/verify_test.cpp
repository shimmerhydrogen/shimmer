/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */
#include "verify_test.h"


bool verify_test(const std::string & name, 
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& vals,
                 const std::vector<double>& ref )
{
    if(vals.size() != ref.size()) 
        return false;

    bool pass = true;

    for (int k = 0; k < vals.size(); ++k)
    {
        // std::cout << vals[k]  <<  std::endl ;
        auto e_val = std::abs((vals[k] - ref.at(k)) /  ref.at(k));
        if(e_val > 1.e-12)
        {
            pass = false;
            break;  
        }
    }
    
    auto passfail = [](bool ipass) {
        return ipass ? "[PASS]" : "[FAIL]";
    };

    std::cout << "  Test " << name << ".........." <<  passfail(pass) << std::endl;

    return pass;
}



bool verify_test(const std::string & name, 
                 const Eigen::SparseMatrix<double>& mat,
                 const std::vector<std::array<double, 3>>& ref )
{
    if(mat.nonZeros() != ref.size()) 
        return false;


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

