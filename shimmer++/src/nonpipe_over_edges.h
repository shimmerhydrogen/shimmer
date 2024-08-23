#pragma once
#include <Eigen/Dense>
#include <memory>
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <tuple>

namespace shimmer
{
namespace edge_station
{
namespace control
{
    enum hardness_type
    {
        HARD,
        SOFT,
        BOUNDARY,
    };

    enum constraint_type
    {
        LOWER,
        GREATER,
        EQUAL,
        LOWER_EQUAL,
        GREATER_EQUAL,
        NONE,
    };

    enum type
    {
        BY_PASS,
        SHUT_OFF,
        POWER_DRIVER,
        PRESSURE_OUT,
        PRESSURE_IN,
        BETA,
        FLUX,
    };

    class constraint
    {
        hardness_type   hardness_;
        constraint_type type_;
        double  value_;

    public:

        constraint(){};
        constraint(const hardness_type& h, const constraint_type& t, double v)
        {};

        bool  check(double variable) const;
        inline double value() const {return value_;};
        inline const constraint_type& type() {return type_;};
    };

    class model
    {
        std::vector<double> coeffs;
        std::vector<int> free_index;
        int control_index_;

        public:
        model() = default;
        model(control::type ctype);
        void set_coefficient(size_t index, double value);
        void set_control_coefficient(double value);
        double coefficient(size_t index);
        double control_coefficient();
    };

    class state
    {
    public:
        model model_;
        control::constraint internal_;
        control::type type_;

        state() = default;
        state(const model& m, const control::constraint& internal):
            model_(m), internal_(internal)
            {};
        state(const control::type& type, const control::constraint& internal):
            model_(model(type)), internal_(internal)
            {};

        virtual bool control_hard(size_t step = 0);

        inline void set_c1(double value){ model_.set_coefficient(0, value);};
        inline void set_c2(double value){ model_.set_coefficient(1, value);};
        inline void set_c3(double value){ model_.set_coefficient(2, value);};
        inline void set_rhs(double value){model_.set_coefficient(3, value);};
        inline double coefficient(size_t index) { return model_.coefficient(index);}
        inline auto control_value() {return internal_.value();};
    };

/*
    class control_power_driver: public state
    {
        double ramp_coeff_;

        public:

        power_driver(double ramp, const constraint& hard):
            type_(control::type::POWER_DRIVER),
            state(type_, hard),
            ramp_coeff_(ramp)
        {};


        bool control_hard(double time)
        {
            // get stored power driver in model
            auto pwd = model_.control_coefficient();

            // 1. Check constraint
            bool pass = internal_.check(pwd);

            // 2. Control variable with hard limit

            // 2.1 Hard limit
            auto pw_nominal = internal_.value();

            // 2.2 Control
            if (!pass)
                model_.set_control_coefficient(pw_nominal);
                // WK: Maybe here I should recheck constraint to be sure new pwd < pwd_nominal;

            return pass;
        }
    };
*/

    auto make_power_driver_control(double PWD_nominal, double ramp);
    auto make_pressure_out_control(double pressure_out_max);
    auto make_pressure_in_control(double pressure_in_min);
    auto make_by_pass_control(const constraint_type& ctype);
    auto make_shutoff_control(const constraint_type& ctype);
    auto make_beta_min_control(double beta_min);
    auto make_beta_max_control(double beta_max);
    auto make_flux_control(double flux_max);

}; //end namespace control

enum external_type
{
    //MIN_GREATER_EQUAL,
    //MAX_LOWER_EQUAL,
    BETA_MIN,
    BETA_MAX,
    P_OUT_MAX,
    P_OUT_MIN,
    P_IN_MIN,
    P_IN_MAX,
    FLUX_MIN,
    FLUX_MAX,
    V_MAX,
    V_MIN,
    PWD_NOMINAL,
    P_THRESHOLD_MIN,
    P_THRESHOLD_MAX,
};


using map_type = std::unordered_map<external_type, control::constraint>;

class station
{

public:
    std::string name_;
    control::state mode_;

    std::vector<control::constraint> internals_;
    std::vector<control::constraint> externals_;

    std::vector<control::state> controls_off;
    std::vector<control::state> controls_on;

    std::vector<bool> active_history_;
    size_t count;

    station() = default;
    station(const std::string& name,
    const std::vector<bool>& active_history,
            const std::vector<control::constraint>& internals,
            const std::vector<control::constraint>& externals):
                    name_(name), active_history_(active_history),
                    internals_(internals), externals_(externals)
        {
        count = 0;
        };

    virtual void activate(size_t step, double target_pressure);

    inline auto which_control_type()
    {
        // size_t idx = count % mode_.size();
        // return mode_[idx].type;
        return control::type::SHUT_OFF;
    };

    bool control_hard(size_t step);
    //bool check_soft(double p, double l, size_t step);

    //set coefficients of the current control model
    void set_c1(double value);
    void set_c2(double value);
    void set_c3(double value);
    void set_rhs(double value);

    inline double model_c1() { return mode_.coefficient(0);};
    inline double model_c2() { return mode_.coefficient(1);};
    inline double model_c3() { return mode_.coefficient(2);};
    inline double model_rhs(){ return mode_.coefficient(3);};

    inline void add_control_on(const control::state& md) { controls_on.push_back(md);};
    inline void add_control_off(const control::state& md){ controls_off.push_back(md);};
};


class compressor : public station
{
    double ramp_coeff_;
    double efficiency;
    int control_node;

    compressor(const std::string& name,
        double efficiency,
        double ramp_coeff,
        const std::vector<bool>& activate_history,
        const std::vector<control::constraint>& internals,
        const std::vector<control::constraint>& externals);

    void activate(size_t step, double target_pressure);
    bool control_hard(size_t step);
    double compute_beta(double pressure_in,
                        double pressure_out);
};


template<typename KEY>
std::unordered_map<KEY,control::constraint>
build_multiple_constraints(const std::vector<std::tuple<KEY, control::constraint_type, double>>& limits,
                            control::hardness_type ht)
{
    std::unordered_map<KEY, control::constraint> constr_vector(limits.size());
    size_t i = 0;
    for(const auto&  l : limits)
        constr_vector[std::get<0>(l)] = control::constraint(ht,
                                                        std::get<1>(l),
                                                        std::get<2>(l));
    return constr_vector;
}


std::vector<control::constraint>
build_multiple_constraints(const std::vector<std::pair<control::constraint_type, double>>& limits,
                            control::hardness_type ht);


auto
make_regulator( const std::vector<bool>& activate_history,
                const std::vector<std::tuple<external_type,
                                             control::constraint_type,
                                             double>> & user_limits);


auto
make_valve( double velocity_limit,
            const std::vector<bool>& activate_history,
            const std::vector<std::tuple<external_type,
                                         control::constraint_type,
                                         double>> & user_limits);


auto
make_compressor(double ramp,
                double efficiency,
                const std::vector<double>& activate_history,
                const std::unordered_map<external_type,
                                        std::pair<control::constraint_type,
                                        double>> & user_limits);


} //end namespace control
}//end namespace shimmer

// 0. Check activation warning
// 1. Check pwd * (step - 1.0) do not get negative
// 2. Use pipes as by_pass....not to heavy? Maybe better just add real stations?
// 3. Marco: In PWD, Ki = ZiTiR.. what about R? at nodes or at pipe? => use c2_nodes, c2_pipes? combination?
// 4. Marco: Maybe here I should recheck constraint to be sure new pwd < pwd_nominal;