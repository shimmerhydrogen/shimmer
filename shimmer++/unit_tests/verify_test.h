/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <Eigen/Sparse>


bool verify_test(const std::string & name, 
                 const Eigen::Matrix<double, Eigen::Dynamic, 1>& vals,
                 const std::vector<double>& ref );


bool verify_test(const std::string & name, 
                 const Eigen::SparseMatrix<double>& mat,
                 const std::vector<std::array<double, 3>>& ref );