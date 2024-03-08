#include "../src/boundary.h"

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
            return (l < (values_[step] + 1.e-14));  
        case L_GREATER_EQUAL:             
            return (l > (values_[step] - 1.e-14));  
        case L_EQUAL:
            return  (std::abs(l - values_[step]) < 1.e-14);
        case P_LOWER_EQUAL:             
            return (p < (values_[step] + 1.e-14));  
        case P_GREATER_EQUAL:             
            return (p > (values_[step] - 1.e-14));  
        case P_EQUAL:
            return (std::abs(p - values_[step]) < 1.e-14);
        default:
            throw std::invalid_argument("Boundary conditions not specified");
    }           
}



state::state(const constraint& b,
             const constraint& i,
             const std::vector<constraint>& e):
             boundary(b), internal(i), externals(e){};



junction::junction()
{
    s0.boundary = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, 0.0);
}


void
remi_wo_backflow::set_boundary(double Pset, const std::vector<pair_input_t>& user_limits_s0)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY,
                                constraint_type::P_EQUAL, Pset);
    auto s0_int = constraint(hardness_type::HARD,
                                constraint_type::L_LOWER_EQUAL, 0.0); 


    std::vector<constraint> s0_ext(user_limits_s0.size());
    size_t i = 0;
    for(const auto&  ec : user_limits_s0) 
        s0_ext[i++] = constraint(hardness_type::SOFT, user_limits_s0[i].first,
                                                    user_limits_s0[i].second);

    states_[0] = state(s0_bnd, s0_int, s0_ext);
}



// WK: This functions are equal since the type of the values changes. Improve this
void
remi_wo_backflow::set_boundary(const vector_t& Pset, const std::vector<pair_input_t>& user_limits_s0)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY,
                                constraint_type::P_EQUAL, Pset);
    auto s0_int = constraint(hardness_type::HARD,
                                constraint_type::L_LOWER_EQUAL, 0.0); 


    std::vector<constraint> s0_ext(user_limits_s0.size());
    size_t i = 0;
    for(const auto&  ec : user_limits_s0) 
        s0_ext[i++] = constraint(hardness_type::SOFT, user_limits_s0[i].first,
                                                    user_limits_s0[i].second);

    states_[0] = state(s0_bnd, s0_int, s0_ext);
}



// WK: This functions are equal since the type of the values changes. Improve this
void 
remi_wo_backflow::set_boundary_to_switch(double Pset,
            const std::vector<pair_input_t>& user_limits_s1)
{
    auto s1_bnd = constraint(hardness_type::BOUNDARY,
                                constraint_type::L_EQUAL, 0.0); 
    auto s1_int = constraint(hardness_type::HARD,
                                constraint_type::P_GREATER_EQUAL, Pset);

    std::vector<constraint> s1_ext(user_limits_s1.size());
    size_t i = 0;
    for(const auto&  ec : user_limits_s1)
        s1_ext[i++] = constraint(hardness_type::SOFT, user_limits_s1[i].first,
                                                    user_limits_s1[i].second);
    states_[1] = state(s1_bnd, s1_int, s1_ext);
}


void 
remi_wo_backflow::set_boundary_to_switch(const vector_t& Pset,
            const std::vector<pair_input_t>& user_limits_s1)
{
    auto s1_bnd = constraint(hardness_type::BOUNDARY,
                                constraint_type::L_EQUAL, 0.0); 

    auto s1_int = constraint(hardness_type::HARD,
                                constraint_type::P_GREATER_EQUAL, Pset);

    std::vector<constraint> s1_ext(user_limits_s1.size());
    size_t i = 0;
    for(const auto&  ec : user_limits_s1)
        s1_ext[i++] = constraint(hardness_type::SOFT, user_limits_s1[i].first,
                                                    user_limits_s1[i].second);
    states_[1] = state(s1_bnd, s1_int, s1_ext);
}


void 
remi_wo_backflow::switch_state()
{
    count_++;
    index_ = count_%num_states_;
    std::cout << "Switch done" << std::endl;
}


bool
remi_wo_backflow::check_hard(double p, double l, size_t step)
{
    if(states_[index_].boundary.check(p, l, step)) 
        return true;

    switch_state();
    return false;
}

bool
remi_wo_backflow::check_soft(double p, double l, size_t step)
{
    bool success = true; 
    for(auto& e : states_[index_].externals)
    {
        if(e.check(p, l, step))
        {
            success = false;
            std::cout << "WARNING: REMI_WO constraint violated => "
                      << e.type() << " " << e.value(step) <<std::endl;                        
        }
    }

    return success;
}

const constraint & remi_wo_backflow::boundary(){ return states_[index_].boundary;}


} //end namespace shimmer