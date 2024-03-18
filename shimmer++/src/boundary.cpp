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
            return (l < (value(step) + 1.e-14));  
        case L_GREATER_EQUAL:             
            return (l > (value(step) - 1.e-14));  
        case L_EQUAL:
            return  (std::abs(l - value(step)) < 1.e-14);
        case P_LOWER_EQUAL:             
            return (p < (value(step) + 1.e-14));  
        case P_GREATER_EQUAL:             
            return (p > (value(step) - 1.e-14));  
        case P_EQUAL:
            return (std::abs(p - value(step)) < 1.e-14);
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
    index_ = count_%num_states_;
}



bool
multiple_states_station::check_hard(double p, double l, size_t step)
{

    if(states_[index_].internal.check(p, l, step)){

        return true;
    } 
    switch_state();

    std::cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<< std::endl;
    std::cout << "WARNING HARD: In "<< name_ <<" constraint violated. SWITCH done." << std::endl;
    std::cout << " * Hard constaint (" << index_ <<") : "<< states_[index_].internal.type() 
                                    << states_[index_].internal.value(step) << std::endl;
    std::cout << " * (press, lrate) : ( " << p <<" , "<< l<< " ) " << std::endl;
    std::cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<< std::endl;

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
    std::cout<< " * index :" << index_ << std::endl;
    std::cout<< " * states.size :" << states_.size() << std::endl;

    return states_[index_].boundary;}



} //end namespace shimmer