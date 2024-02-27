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
    virtual void set_boundary(const vector_t& vals) = 0;
    virtual void set_hard_constraints(const std::vector<pair_input_t>&) = 0;
    virtual void set_soft_constraints(const std::vector<pair_input_t>&) = 0;

    virtual bool check_hard(double, double, size_t) = 0;
    virtual bool check_soft(double, double, size_t) = 0;
    virtual const constraint& boundary() = 0;
    virtual void print(){std::cout << "STATION" << std::endl;}
};



class junction: public station
{
    state s0;

public:
    junction();

    inline void set_boundary(const vector_t& vals){}
    inline void set_hard_constraints(const std::vector<pair_input_t>&){};
    inline void set_soft_constraints(const std::vector<pair_input_t>&){};
    inline bool check_hard(double, double, size_t) {return true;};
    inline bool check_soft(double, double, size_t) {return true;};
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "JUNCTION" << std::endl;}
};



class inlet_station: public station
{
    state s0;
public:
    inline inlet_station(){};
    inline void set_boundary(const vector_t& vals)
    {
        s0.boundary = constraint(hardness_type::BOUNDARY,
                                 constraint_type::P_EQUAL, vals); 
    }
    inline void set_hard_constraints(const std::vector<pair_input_t>&){};
    inline void set_soft_constraints(const std::vector<pair_input_t>&){};
    inline bool check_hard(double p, double l, size_t step){return true;};
    inline bool check_soft(double p, double l, size_t step){return true;};    
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "INLET" << std::endl;}
};


class outlet_station: public station
{
    state s0;

public:
    inline outlet_station(){};   
    inline void set_boundary(const vector_t& vals)
    {
        s0.boundary = constraint(hardness_type::BOUNDARY,
                                 constraint_type::L_EQUAL, vals); 

    }  

    inline void set_hard_constraints(const std::vector<pair_input_t>&){};
    inline void set_soft_constraints(const std::vector<pair_input_t>&){};
    inline bool check_hard(double p, double l, size_t step){return true;};
    inline bool check_soft(double p, double l, size_t step){return true;};    
    inline const constraint& boundary(){return s0.boundary;};
    inline void print(){std::cout << "OUTLET" << std::endl;}
};


/*
class remi_wo_backflow: public station
{
    size_t count_;
    size_t index_;
    size_t num_states_;
    std::vector<state> states_;

public:
    remi_wo_backflow(){};
    remi_wo_backflow(double ,
                     const std::vector<pair_input_t>& ,
                     const std::vector<pair_input_t>& );
    remi_wo_backflow(const vector_t& ,
                     const std::vector<pair_input_t>& ,
                     const std::vector<pair_input_t>& );

    void switch_state();
    bool check_hard(double, double, size_t);
    bool check_soft(double, double, size_t);
    const constraint & boundary();
};
*/



} //end namespace shimmer