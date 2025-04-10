#pragma once
#include <Eigen/Dense>
#include <memory>
#include <iostream>
#include <vector>
#include <cassert>
namespace shimmer {

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

enum hardness_type
{
    HARD,
    SOFT,
    BOUNDARY,
};


enum constraint_type
{   
    L_LOWER_EQUAL,
    L_GREATER_EQUAL,
    L_EQUAL,
    P_LOWER_EQUAL,
    P_GREATER_EQUAL,
    P_EQUAL,
};
std::ostream& operator <<(std::ostream& os, const constraint_type& o);


using pair_input_t = std::pair<constraint_type,double>;


class constraint
{
    hardness_type   hardness_;
    constraint_type type_;

    bool is_varying_in_time; 
    double  value_;
    vector_t values_;

public:

    inline  constraint(){};
    constraint(const hardness_type& h, const constraint_type& t, double v);
    constraint(const hardness_type& h, const constraint_type& t, const vector_t& v);

    bool  check(double p, const double l, size_t step);
    inline double value() const {
        assert(is_varying_in_time && "Values are varying in time and need to provide a step.");
        return value_;
    }
    inline double value(const size_t& step) const 
    {   
        if(is_varying_in_time)
            return values_[step];
        return value_;
    };
    inline const constraint_type& type() {return type_;};
};



struct state
{
    constraint  boundary;  
    constraint  internal;
    std::vector<constraint> externals;

    inline state(){};
    inline state(const constraint& b,const constraint&i, const std::vector<constraint>& e):
             boundary(b), internal(i), externals(e){};
};



class station
{
public: 
    std::string name_;
    station(){};
    virtual ~station() {}

    virtual void set_state(const state& s) = 0;   
    virtual void set_state_to_switch(const state& s){};

    virtual bool check_hard(double, double, size_t) {return true;};
    virtual bool check_soft(double, double, size_t) {return true;};
    virtual const constraint& boundary() = 0;
    virtual void print(){std::cout << "STATION" << std::endl;}
    virtual void switch_state(){};
};

std::vector<constraint>
build_user_constraints(const std::vector<pair_input_t>& user_limits);


class junction: public station
{
    state s0;
public:
    junction()
    { 
        s0.boundary = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, 0.0);
    }

    inline void set_state(const state& s){};
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "JUNCTION" << std::endl;}
};



class one_state_station: public station
{
    state s0;

public:
    one_state_station(std::string name) {name_ = name;}
    inline void set_state(const state& s){s0 = s;}   
 
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << name_ << std::endl;}
};

/* Make exit L-regulated station */
template<typename VALUE>
one_state_station
make_station_exit_l_reg(const VALUE& vals, const std::vector<pair_input_t>& user_limits)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, vals); 
    auto s0_int = constraint(hardness_type::HARD, constraint_type::L_GREATER_EQUAL, 0.0); 
    auto s0_ext = build_user_constraints(user_limits);

    auto s0 = state(s0_bnd, s0_int, s0_ext);

    one_state_station consumption_station("CONSUMPTION_WO_PRESSURE"); 
    consumption_station.set_state(s0);
    
    return consumption_station;
}


namespace priv {

template<typename VALUE>
one_state_station
make_station_inlet(const VALUE& vals)
{
    state s0;
    s0.boundary = constraint(hardness_type::BOUNDARY, constraint_type::P_EQUAL, vals);

    one_state_station inlet("ENTRY");
    inlet.set_state(s0);
    return inlet;
}

template<typename VALUE>
one_state_station
make_station_outlet(const VALUE& vals)
{
    state s0;
    s0.boundary = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, vals); 

    one_state_station out("EXIT"); 
    out.set_state(s0);
    return out;
}

} // namespace priv

class multiple_states_station: public station
{
    bool only_one_switch_; 
    size_t count_;
    size_t index_;
    size_t num_states_;
    std::vector<state> states_;

public:
    multiple_states_station(std::string name,
                            bool one_switch_allowed = false):
    only_one_switch_(one_switch_allowed)
    {
        name_ = name;
        count_ = 0;
        index_ = 0;
        num_states_ = 2;
    };

    multiple_states_station(std::string name, size_t num_states,
                            bool one_switch_allowed = false):
    only_one_switch_(one_switch_allowed),
    num_states_(num_states)
    {
        name_ = name;
        count_ = 0;
        index_ = 0;
    };


    inline void set_state(const state& s)
    {
        assert((states_.size() < num_states_)&& "Station have more states than given by definition");
        states_.push_back(s);
    }
    void set_state_to_switch(){};
    void switch_state();
    bool check_hard(double, double, size_t);
    bool check_soft(double, double, size_t);
    const constraint & boundary();
    inline void print(){std::cout << name_  << std::endl;}
};


/* Make entry L-regulated station */
template<typename L_TYPE,typename P_TYPE >
multiple_states_station
make_station_entry_l_reg(double factor, const P_TYPE& Pset, const L_TYPE& Lset,
                         const std::vector<pair_input_t>& user_limits_s0,
                         const std::vector<pair_input_t>& user_limits_s1,
                    bool one_switch_allowed = false)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, Lset);
    auto s0_int = constraint(hardness_type::HARD,constraint_type::P_LOWER_EQUAL, factor * Pset); 
    auto s0_ext = build_user_constraints(user_limits_s1);

    auto s1_bnd = constraint(hardness_type::BOUNDARY, constraint_type::P_EQUAL, Pset); 
    auto s1_int = constraint(hardness_type::HARD, constraint_type::L_LOWER_EQUAL, 0.0);
    auto s1_ext = build_user_constraints(user_limits_s1);

    auto s0 = state(s0_bnd, s0_int, s0_ext);
    auto s1 = state(s1_bnd, s1_int, s1_ext);

    multiple_states_station inj_station("INJECTION_WP_CONTROL", one_switch_allowed);

    inj_station.set_state(s0);
    inj_station.set_state(s1);

    return inj_station;
}

/* Make entry P-regulated station */
template<typename VALUE_TYPE>
multiple_states_station
make_station_entry_p_reg(  const VALUE_TYPE& Pset_s0,
                           const std::vector<pair_input_t>& user_limits_s0,
                           const std::vector<pair_input_t>& user_limits_s1,
                           bool one_switch_allowed = false)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY, constraint_type::P_EQUAL, Pset_s0);
    auto s0_int = constraint(hardness_type::HARD,constraint_type::L_LOWER_EQUAL, 0.0); 
    auto s0_ext = build_user_constraints(user_limits_s0);

    // THIS MUST BE MODFIED SINCE IT WAS SET TO ZERO ONLY DUE TO COMPARISON WITH MATLAB CODE
    auto s1_bnd = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, 0.0); //Carefull here!!!!! 
    auto s1_int = constraint(hardness_type::HARD, constraint_type::P_GREATER_EQUAL, 0.0);
    auto s1_ext = build_user_constraints(user_limits_s1);

    auto s0 = state(s0_bnd, s0_int, s0_ext);
    auto s1 = state(s1_bnd, s1_int, s1_ext);

    multiple_states_station remi("REMI_WO_BACKFLOW", one_switch_allowed);
    remi.set_state(s0);
    remi.set_state(s1);

    return remi;
}

} //end namespace shimmer