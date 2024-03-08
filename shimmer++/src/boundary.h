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
    state(const constraint&,const constraint&, const std::vector<constraint>& e);
};


class station
{
public: 

    station(){};
    virtual void set_boundary(const vector_t& vals, const std::vector<pair_input_t>& user_limits ={}) = 0;
    virtual void set_boundary(double val, const std::vector<pair_input_t>& user_limits = {} ) = 0;
    virtual void set_boundary_to_switch(const vector_t& vals, const std::vector<pair_input_t>& user_limits = {}){};
    virtual void set_boundary_to_switch(double val, const std::vector<pair_input_t>& user_limits = {} ){};

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
    junction();

    inline void set_boundary(const vector_t& vals,const std::vector<pair_input_t>& user_limits ={}){}
    inline void set_boundary(double val,  const std::vector<pair_input_t>& user_limits ={}){}
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "JUNCTION" << std::endl;}
};



class inlet_station: public station
{
    state s0;
public:
    inline inlet_station(){};
    inline void set_boundary(const vector_t& vals,const std::vector<pair_input_t>& user_limits ={})
    {
        s0.boundary = constraint(hardness_type::BOUNDARY,
                                 constraint_type::P_EQUAL, vals); 
    }
    inline void set_boundary(double val,const std::vector<pair_input_t>& user_limits ={})    
    {
        s0.boundary = constraint(hardness_type::BOUNDARY,
                                 constraint_type::P_EQUAL, val); 
    }

    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "INLET" << std::endl;}
};


class outlet_station: public station
{
    state s0;

public:
    inline outlet_station(){};   
    inline void set_boundary(const vector_t& vals, const std::vector<pair_input_t>& user_limits ={})
    {
        s0.boundary = constraint(hardness_type::BOUNDARY,
                                 constraint_type::L_EQUAL, vals); 

    }  
    inline void set_boundary(double val,const std::vector<pair_input_t>& user_limits ={})
    {
        s0.boundary = constraint(hardness_type::BOUNDARY,
                                 constraint_type::L_EQUAL, val); 

    }  
 
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "OUTLET" << std::endl;}
};



class remi_wo_backflow: public station
{
    size_t count_;
    size_t index_;
    size_t num_states_;
    std::vector<state> states_;

public:
    remi_wo_backflow(){
        count_ = 0;
        num_states_ = 2;
        states_.resize(num_states_);
    };

    void set_boundary(double,
                     const std::vector<pair_input_t>&);
    void set_boundary(const vector_t& ,
                     const std::vector<pair_input_t>&);
    void set_boundary_to_switch(double,
                     const std::vector<pair_input_t>&);
    void set_boundary_to_switch(const vector_t& ,
                     const std::vector<pair_input_t>&);

    void switch_state();
    bool check_hard(double, double, size_t);
    bool check_soft(double, double, size_t);
    const constraint & boundary();
    inline void print(){std::cout << "REMI WO BACKFLOW" << std::endl;}
};




} //end namespace shimmer