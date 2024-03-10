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


using pair_input_t = std::pair<constraint_type,double>;


std::ostream& operator <<(std::ostream& os, const constraint_type& o);


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

    station(){};

    virtual void set_state(const state& s) = 0;   
    virtual void set_state_to_switch(const state& s){};   

    virtual bool check_hard(double, double, size_t) {return true;};
    virtual bool check_soft(double, double, size_t) {return true;};
    virtual const constraint& boundary() = 0;
    virtual void print(){std::cout << "STATION" << std::endl;}
    virtual void switch_state(){};
};


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



template<typename VALUE>
state
make_inlet(const VALUE& vals)
{
    state s0;
    s0.boundary = constraint(hardness_type::BOUNDARY, constraint_type::P_EQUAL, vals);
    return s0;
}


class inlet_station: public station
{
    state s0;
public:
    inline inlet_station(){};
    inline void set_state(const state& s){s0 = s;};
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "INLET" << std::endl;}
};


std::vector<constraint>
build_user_constraints(const std::vector<pair_input_t>& user_limits);


template<typename VALUE>
state
make_consumption_wo_press(const VALUE& Lset, const std::vector<pair_input_t>& user_limits_s0)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, Lset); 
    auto s0_int = constraint(hardness_type::HARD, constraint_type::L_GREATER_EQUAL, 0.0); 
    auto s0_ext = build_user_constraints(user_limits_s0);

    return state(s0_bnd, s0_int, s0_ext);
}



template<typename VALUE>
state
//make_consumption_wo_press(const vector_t& vals)
make_outlet(const VALUE& vals)
{
    state s0;
    s0.boundary = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, vals); 
    return s0;
}



class outlet_station: public station
{
    state s0;

public:
    inline outlet_station(){};
    inline void set_state(const state& s){s0 = s;}   
 
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "OUTLET" << std::endl;}
};



template<typename VALUE_TYPE>
std::pair<state, state>
make_remi_wo_backflow(  const VALUE_TYPE& Pset,
                        const std::vector<pair_input_t>& user_limits_s0,
                        const std::vector<pair_input_t>& user_limits_s1)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY, constraint_type::P_EQUAL, Pset);
    auto s0_int = constraint(hardness_type::HARD,constraint_type::L_LOWER_EQUAL, 0.0); 
    auto s0_ext = build_user_constraints(user_limits_s0);

    auto s1_bnd = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, 0.0); 
    auto s1_int = constraint(hardness_type::HARD, constraint_type::P_GREATER_EQUAL, Pset);
    auto s1_ext = build_user_constraints(user_limits_s1);

    auto s0 = state(s0_bnd, s0_int, s0_ext);
    auto s1 = state(s1_bnd, s1_int, s1_ext);

    return std::make_pair(s0, s1);
}



class remi_wo_backflow: public station
{
    size_t count_;
    size_t index_;
    size_t num_states_;
    std::vector<state> states_;

public:
    remi_wo_backflow(){
        count_ = 0;
        index_ = 0;
        num_states_ = 2;
        states_.resize(num_states_);
    };

    inline void set_state(const state& s){states_[0] = s;}
    inline void set_state_to_switch(const state& s){states_[1] = s;};
    void switch_state();
    bool check_hard(double, double, size_t);
    bool check_soft(double, double, size_t);
    const constraint & boundary();
    inline void print(){std::cout << "REMI WO BACKFLOW" << std::endl;}
};



template<typename L_TYPE,typename P_TYPE >
std::pair<state, state>
make_inj_w_pressure(double factor, const P_TYPE& Pset, const L_TYPE& Lset,
                    const std::vector<pair_input_t>& user_limits_s0,
                    const std::vector<pair_input_t>& user_limits_s1)
{
    auto s0_bnd = constraint(hardness_type::BOUNDARY, constraint_type::L_EQUAL, Lset);
    auto s0_int = constraint(hardness_type::HARD,constraint_type::P_LOWER_EQUAL, factor * Pset); 
    auto s0_ext = build_user_constraints(user_limits_s1);

    auto s1_bnd = constraint(hardness_type::BOUNDARY, constraint_type::P_EQUAL, Pset); 
    auto s1_int = constraint(hardness_type::HARD, constraint_type::L_LOWER_EQUAL, 0.0);
    auto s1_ext = build_user_constraints(user_limits_s1);

    auto s0 = state(s0_bnd, s0_int, s0_ext);
    auto s1 = state(s1_bnd, s1_int, s1_ext);

    return std::make_pair(s0,s1);
}



class injection_wp_control: public station
{
    size_t count_;
    size_t index_;
    size_t num_states_;
    std::vector<state> states_;

public:
    injection_wp_control(){
        count_ = 0;
        index_ = 0;
        num_states_ = 2;
        states_.resize(num_states_);
    };

    inline void set_state(const state& s){states_[0] = s;}
    inline void set_state_to_switch(const state& s){states_[1] = s;};
    void switch_state();
    bool check_hard(double, double, size_t);
    bool check_soft(double, double, size_t);
    const constraint & boundary();
    inline void print(){std::cout << "INJECTION W PRESSURE" << std::endl;}
};




} //end namespace shimmer