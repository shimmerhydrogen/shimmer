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

#include "../src/boundary.h"
#include <fstream>
namespace shimmer {

std::ostream& operator <<(std::ostream& os, const constraint_type& o)
{
    switch(o)
    {
        case (constraint_type::L_LOWER_EQUAL):
            os << "Flux rate <=";
            return os;
        case (constraint_type::L_GREATER_EQUAL):
            os << "Flux rate >=";
            return os;
        case (constraint_type::L_EQUAL):
            os << "Flux rate =";
            return os;
        case (constraint_type::P_LOWER_EQUAL):
            os << "Pressure <=";
            return os;
        case (constraint_type::P_GREATER_EQUAL):
            os << "Pressure >=";
            return os;
        case (constraint_type::P_EQUAL):
            os << "Pressure =";
            return os;
        default:
            throw std::invalid_argument("ERROR: constraint type not found.");
    }
}



constraint::constraint(const hardness_type& h, const constraint_type& t, double v):
            hardness_(h), type_(t), value_(v)
{
    if(h == hardness_type::BOUNDARY)
    {
        if (t != constraint_type::L_EQUAL && t != constraint_type::P_EQUAL)
            throw std::invalid_argument("Boundary built with incorrect constraint.");
    }

    is_varying_in_time = false;
}



constraint::constraint(const hardness_type& h, const constraint_type& t, const vector_t& v):
            hardness_(h), type_(t), values_(v)
{
    if(h == hardness_type::BOUNDARY)
    {
        if (t != constraint_type::L_EQUAL && t != constraint_type::P_EQUAL)
            throw std::invalid_argument("Boundary built with incorrect constraint.");
    }
    is_varying_in_time = true;
}



bool
constraint::check(double p, double l, size_t step)
{
    switch(type_)
    {
        case L_LOWER_EQUAL:   
            return (l <= value(step));  
        case L_GREATER_EQUAL: 
            return (l >= value(step) );  
        case L_EQUAL:
            return ((std::abs(l - value(step))/value(step)) < 1.e-15);
        case P_LOWER_EQUAL:             
            return (p < value(step));  
        case P_GREATER_EQUAL:             
            return (p >= value(step));  
        case P_EQUAL:
            return ((std::abs(p - value(step))/value(step)) < 1.e-15);
        default:
            throw std::invalid_argument("Boundary conditions not specified");
    }           
}



std::vector<constraint>
build_user_constraints(const std::vector<pair_input_t>& user_limits)
{
   std::vector<constraint> externals(user_limits.size());
    size_t i = 0;
    for(const auto&  ec : user_limits) 
        externals[i++] = constraint(hardness_type::SOFT, 
                                user_limits[i].first,
                                user_limits[i].second);
    return externals;
}



void 
multiple_states_station::switch_state()
{
    count_++;
    
    /*  
        Only_one_switch_ is added only for testing reasons against 
        the Matlab code. Must be set as FALSE by default
    */
    index_ = (only_one_switch_)? 1 :  
                                 count_%num_states_;
}



bool
multiple_states_station::check_hard(double p, double l, size_t step)
{

    if(states_[index_].internal.check(p, l, step))
        return true;

    std::cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<< std::endl;
    std::cout << "WARNING HARD: In "<< name_ <<" constraint violated. SWITCH done." << std::endl;
    std::cout << " * Hard constaint (" << index_ <<") : "<< states_[index_].internal.type() 
                                    << states_[index_].internal.value(step) << std::endl;
    std::cout << " * (press, lrate) : ( " << p <<" , "<< l<< " ) " << std::endl;
    std::cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<< std::endl;

    switch_state();

    return false;
}



bool
multiple_states_station::check_soft(double p, double l, size_t step)
{
    bool success = true; 
    for(auto& e : states_[index_].externals)
    {
        if(!e.check(p, l, step))
        {
            success = false;
            std::cout << "WARNING SOFT:" << name_  << " constraint violated." << std::endl;
            std::cout << " * Soft constraint ("<<index_ <<") : " << e.type() << " " << e.value(step) <<std::endl;     
            std::cout << " * (press , lrate) :  ("<<p << "," << l<< ") " << std::endl;                        
        }
    }

    return success;
}


const constraint & multiple_states_station::boundary(){ 
    //std::cout<< " * index :" << index_ << std::endl;
    //std::cout<< " * states.size :" << states_.size() << std::endl;

    return states_[index_].boundary;}



} //end namespace shimmer