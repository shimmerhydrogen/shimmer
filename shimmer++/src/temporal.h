
/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * Karol Cascavita (C) 2023
 * karol.cascavita@polito.it  
 */

#pragma once

#include <Eigen/Sparse>

namespace shimmer{

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>; 


struct variable
{
    vector_t pressure;
    vector_t flux;
    vector_t L_rate;

    variable();
    variable(const vector_t&p,const vector_t&f,const vector_t&l);
    vector_t make_vector();
};

}
