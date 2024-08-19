#pragma once
#include <Eigen/Dense>
#include <memory>
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>

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

        inline  constraint(){};
        constraint(const hardness_type& h, const constraint_type& t, double v)
        {};

        bool  check(double variable, size_t step)
        {
            switch (type_)
            {
                case LOWER_EQUAL:
                    return (variable < (value_ + 1.e-14));
                case GREATER_EQUAL:
                    return (variable > (value_ - 1.e-14));
                case EQUAL:
                    return (std::abs(variable - value_) < 1.e-14);
                case LOWER:
                    return (variable < value_);
                case GREATER:
                    return (variable > value_);
                case NONE:
                    return true;
                default:
                    throw std::invalid_argument("Boundary conditions not specified");
            }
        };

        inline double value() const {return value_;}

        inline const constraint_type& type() {return type_;};
    };

    class model
    {
        std::vector<double> coeffs;
        std:vector<int> free_index;
        int control_index_;

        public:
        model(control::type ctype)
        {
            switch (ctype)
            {
            case control::type::BY_PASS:
                coeffs[0] =-1.0;
                coeffs[1] = 1.0;
                coeffs[2] = 0.0;
                coeffs[3] = 0.0;
                control_index_ = -1;
                break;
            case control::type::SHUT_OFF:
                coeffs[0] = 0.0;
                coeffs[1] = 0.0;
                coeffs[2] = 1.0;
                coeffs[3] = 0.0;
                control_index_ = -1;
                break;
            case control::type::POWER_DRIVER:
                free_index = std : vector<int>{ 0,1,2,3 };
                control_index_ = 3;
                break;
            case control::type::PRESSURE_IN:
                coeffs[0] = 1.0;
                coeffs[1] = 0.0;
                coeffs[2] = 0.0;
                free_index = std:vector<int>{ 3 };
                control_index_ = 3;
                break;
            case control::type::PRESSURE_OUT:
                coeffs[0] = 0.0;
                coeffs[1] = 1.0;
                coeffs[2] = 0.0;
                free_index = std:vector<int>{ 3 };
                control_index_ = 3;
                break;
            case control::type::FLUX:
                coeffs[0] = 0.0;
                coeffs[1] = 0.0;
                coeffs[2] = 1.0;
                free_index = std:vector<int>{ 3 };
                control_index_ = 3;
                break;
            case control::type::BETA:
                coeffs[1] = 1.0;
                coeffs[2] = 0.0;
                coeffs[3] = 0.0;
                free_index = std:vector<int>{0};
                control_index_ = 1;
                break;
            default:
                break;
            }
        }

        void set_coefficient(size_t index, double value)
        {
            assert(free_index.size() > 0);
            assert("Error: In model of edge station, coefficient at index is not modifiable" && std::find(free_index.begin(), free_index.end(), item) != vec.end())

            coeffs[index] = value;
        }
        void set_control_coefficient(double value)
        {
            coeffs[control_index_] = value;
        }

        double control_coefficient()
        {
            if (control_index_ < 0)
                return INF;

            return coeffs[control_index_];
        }
    };

    class state
    {
    public:
        model model_;
        control::constraint internal_;
        control::type type_;

        state(const model& m, const control::constraint& internal):
            model_(m), internal_(internal)
            {};
        state(const control::type& type, const control::constraint& internal):
            model_(model(type)), internal_(internal)
            {};

        virtual bool control_hard(double time = 0)
        {
            // Get stored variable in model
            auto variable = model_.control_coefficient();

            // 1. Check constraint
            bool pass = internal_.check(variable);

            // 2. Control variable with hard limit
            if (!pass)
                model_.set_control_coefficient(internal_.value());

            return pass;
        }

        inline void set_c1(double value){ m.set_cofficient(0, value)};
        inline void set_c2(double value){ m.set_cofficient(1, value)};
        inline void set_c3(double value){ m.set_cofficient(2, value)};
        inline void set_rhs(double value)
        {
            m.set_cofficient(3, value)
        };

        inline auto control_value() {return internal_.value()};
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
            // Get stored power driver in model
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

    auto
    make_power_driver_control(double PWD_nominal, double ramp)
    {
        auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                           constraint_type::GREATER_EQUAL,
                                           PWD_nominal);
        return state(control::type::POWER_DRIVER, internal);
        //return control_power_driver(ramp, internal);
    }

    auto
    make_pressure_out_control(double pressure_out_max)
    {
        auto internal = control::constraint(hardness_type::HARD_CONTROL,
                    constraint_type::LOWER_EQUAL, pressure_out_max);

        return state(control::type::PRESSURE_OUT, internal);
    }

    auto
    make_pressure_in_control(double pressure_in_min)
    {
        auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                           constraint_type::GREATER_EQUAL,
                                           pressure_in_min);

        return state(control::type::PRESSURE_IN, internal);
    }

    auto
    make_by_pass_control(const constraint_type& ctype)
    {
        auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                        ctype, 0);

        return state(control::type::BY_PASS, internal);
    }

    auto
    make_shutoff_control(const constraint_type& ctype)
    {
        auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                        ctype, 1);

        return  state(control::type::SHUT_OFF, internal);
    }

    auto
    make_beta_min_control(double beta_min)
    {
        auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                        constraint_type::LOWER_EQUAL, beta_min);

        return  state(control::type::BETA, internal);
    }

    auto
    make_beta_max_control(double beta_max)
    {
        auto internals = control::constraint(hardness_type::HARD_CONTROL,
                                        constraint_type::GREATER_EQUAL, beta_max);

        return  state(control::type::BETA, internal);
    }

    auto
    make_flux_control(double flux_max)
    {
        auto internals = control::constraint(hardness_type::HARD_CONTROL,
                                        constraint_type::GREATER_EQUAL, flux_max);

        return  state(control::type::FLUX, internal);
    }

} //end namespace control

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
        FLUX_MINL,
        FLUX_MAX,
        V_MAX,
        V_MIN,
        PWD_NOMINAL,
        P_THRESHOLD_MIN,
        P_THRESHOLD_MAX,
    }


    class station
    {
        size_t count;
        std::string name_;
        bool is_on_;
        std::vector<bool> active_history_;

        std::vector<control::constraint> internals_;
        std::unordered_map<control::constraint> externals_;

        std::vector<control::state> controls_off;
        std::vector<control::state> controls_on;

        station(const std::string& name,
                const std::vector<bool>& active_history,
                const std::vector<control::constraint>& internals,
                const std::unordered_map<control::constraint>& externals):
                     name_(name), active_history_(active_history),
                     internals_(internals), externals_(externals)
        {
            count = 0;
        };

        inline virtual auto is_active(size_t step, const infrastructure_graph& graph, const variable& var)
            {return active_history_[step]};

        inline auto which_control_type()
        {
            size_t idx = count % controls.size();
            return controls[idx].type;
        };

        bool
        control_hard(double time)
        {
            size_t idx = count % controls.size();
            bool pass = controls[idx].control_hard(time);

            if (!pass)
            {
                count++;
                return false;
            }
        }

            /*
                bool
        control_hard(double time)
        {
            size_t idx = count%controls.size();

            for (size_t i = 0; i < controls.size(); i++; idx++)
            {
                bool pass = controls[idx].control_hard(time);

                if (!pass)
                {
                    count++;
                    return false;
                }
            }

            return true;
        }
        */

        bool
        check_soft(double p, double l, size_t step)
        {
            bool success = true;
            for(auto& e : externals_)
            {
                if(!e.check(vars, step))
                {
                    success = false;
                    std::cout << "WARNING SOFT:" << name_  << " constraint violated." << std::endl;
                    std::cout << " * Soft constraint (-) : " << e.type() << " " << e.value(step) <<std::endl;
                    std::cout << " * (press , lrate) :  ("<<p << "," << l<< ") " << std::endl;
                }
            }

            return success;
        }

        //set coefficients of the current control model
        void set_c1(double value)
        {
            size_t idx = count % controls.size();
            controls[idx].set_coefficient(0, value);
        }
        void set_c2(double value)
        {
            size_t idx = count % controls.size();
            controls[idx].set_coefficient(1, value);
        }
        void set_c3(double value)
        {
            size_t idx = count % controls.size();
            controls[idx].set_coefficient(2, value);
        }
        void set_rhs(double value)
        {
            size_t idx = count % controls.size();
            controls[idx].set_coefficient(3, value);
        };
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
                    const std::unordered_map<control::constraint>& externals):
            station(name, activate_history, internals, externals){};

        bool control_hard(size_t step, const infrastructure_graph& graph, const variable& var)
        {
            bool is_on = is_active(step, graph, var);

            auto controls = (is_on) ? std::make_shared<std::vector<control::state>>(controls_on)
                                    : std::make_shared<std::vector<control::state>>(controls_off);

            size_t idx = count % controls.size();
            bool pass = controls[idx].control_hard(time);

            if (!pass)
            {
                count++;
                return false;
            }

            if(controls[idx].type() == control::type::POWER_DRIVER)
            {
                auto pwd = controls[idx].model_.control_coefficient();
                auto pwd_nominal = controls[idx].internal_value();
                auto pwd_control = pwd * (step - 1.0) + ramp_coeff_ * pwd_nominal;
                model_.set_control_coefficient(pwd_control);
            }

            return true;
        }

        bool is_active(size_t step, const infrastructure_graph& graph, const variable& var)
        {
            /* Warning: control node is set for the moment with itor = begin(),
                        but it must be updated to source. This can be done easily
                        with the new infrastructure already in branch MC/db_interface
                        (a1847d9) and update to source once all is integrated
            */
            auto itor = boost::edges(g).first;
            auto s = source(*itor, g);
            auto control_node = g[s].node_num;

            auto pn = var_.pressure[control_node];
            bool pass_down = externals_[BETA_MIN].check(pn);
            bool pass_up   = externals_[BETA_MAX].check(pn);

            bool is_on   = (!active_history_[step] & pass_up)
                             || (active_history_[step] & !pass_down);

            return is_on;
        }

    };


    template<typename CONSTR_NONPIPE>
    std::unordered_map<CONSTR_NONPIPE>
    build_multiple_constraints(const std::vector<std::tuple<external_type, constraint_type, double>>& limits,
                                hardness_type ht)
    {
        std::unordered_map<CONSTR_NONPIPE> constr_vector(limits.size());
        size_t i = 0;
        for(const auto&  l : limits)
            constr_vector[get<0>(limits)] = constraint(ht,
                                            get<1>(limits),
                                            get<2>(limits));
        return constr_vector;
    }

    template<typename CONSTR_NONPIPE>
    std::vector<CONSTR_NONPIPE>
    build_multiple_constraints(const std::vector<std::pair<constraint_type, double>>& limits,
                                hardness_type ht)
    {
        std::vector<CONSTR_NONPIPE> constr_vector(limits.size());
        size_t i = 0;
        for(const auto&  l : limits)
            constr_vector[i++] = constraint(ht,
                                            limits[i].first,
                                            limits[i].second);
        return constr_vector;
    }

    auto
    make_regulator(const std::vector<bool>& activate_history,
                   const std::vector<std::tuple<external_type, constraint_type, double>>& user_limits)
    {
        auto internal = control::constraint(hardness_type::HARD,
                                        constraint_type::GREATER_EQUAL,
                                         velocity_limit);
        auto externals = build_multiple_constraints<control::constraint>(user_limits);

        station regulator("REGULATOR", activate_history, internal, externals);

        auto c_by_pass = make_by_pass_control(constraint_type::GREATER_EQUAL);
        auto c_shutoff = make_shutoff_control(constraint_type::LOWER);

        valve.set_control_on(c_by_pass);
        valve.set_control_off(c_shutoff);

        return;
    }

    auto
    make_valve(const std::vector<bool>& activate_history,
               double velocity_limit,
               const std::vector<std::tuple<external_type, constraint_type, double>>& user_limits)
    {
        auto internal = control::constraint(hardness_type::HARD,
                                        constraint_type::GREATER_EQUAL,
                                         velocity_limit);
        auto externals = build_multiple_constraints<control::constraint>(user_limits);

        station valve("VALVE", activate_history, internal, externals);

        auto c_by_pass = make_by_pass_control(constraint_type::NONE);
        auto c_shutoff = make_shutoff_control(constraint_type::NONE);

        valve.set_control_on(c_by_pass);
        valve.set_control_off(c_shutoff);

        return;
    }

    auto
    make_compressor(size_t control_node,
                    double ramp,
                    double efficiency,
                    double power_driver_nominal,
                    const std::vector<double>& activate_history,
                    const std::vector<std::pair<constraint_type, double>>& hard_limits,
                    const std::vector<std::tuple<external_type, constraint_type, double>>& user_limits)
    {
        auto thresholds = build_multiple_constraints<constraint>(pressure_limits, hardness_type::HARD);
        auto internals  = build_multiple_constraints<control::constraint>(hard_limits, hardness_type::HARD);
        auto externals  = build_multiple_constraints<control::constraint>(user_limits, hardness_type::SOFT);

        compressor cmp("COMPRESSOR", activate_history, internals, externals, thresholds);

        cmp.set_control_node(control_node, thresholds);

        auto c_by_pass = make_by_pass_control(constraint_type::EQUAL);
        auto c_shutoff = make_shutoff_control(constraint_type::EQUAL);

        cmp.set_control_off(c_by_pass);
        cmp.set_control_off(c_shutoff);

        const auto beta_min  = externals[BETA_MIN].value();
        const auto beta_max  = externals[BETA_MAX].value();
        const auto p_in_min  = externals[P_IN_MIN].value();
        const auto p_out_max = externals[P_OUT_MAX].value();
        const auto flux_max  = externals[FLUX_MAX].value();

        auto c_power_driver = make_power_driver_control(power_driver_nominal, ramp);
        auto c_beta_min = make_beta_min_control(beta_min);
        auto c_beta_max = make_beta_max_control(beta_max);
        auto c_press_in = make_pressure_in_control(pressure_in_min);
        auto c_press_out = make_pressure_out_control(pressure_out_max);
        auto c_flux = make_flux_control(flux_max);

        cmp.set_control_on(c_power_driver);
        cmp.set_control_on(c_press_in);
        cmp.set_control_on(c_press_out);
        cmp.set_control_on(c_beta_max);
        cmp.set_control_on(c_beta_min);
        cmp.set_control_on(c_flux);

        return cmp;
    }

/*
    double
    compressor_power_driver(double beta,
                            double pressure_in,
                            double pressure_out,
                            double G,
                            double efficiency,
                            double ZTR_inlet,
                            double adiabatic_index = 1.4)
    {
        // 4. Set beta
        auto beta = set_beta(p_in, p_out, internals);
        auto K = ZTR_inlet/efficiency;
        auto ck = (adiabatic_index - 1.0) / adiabatic_index;

        // 1. Compute power driver
        return (K * G / ck) * (std::pow(beta, ck) - 1.0);
    }
*/

    double
    compressor_beta( double p_in,
                     double p_out,
                     const std::unordered_map<control::constraint>& internals)
    {
        double press_rate = p_out / p_in;


        const auto& beta_min_constr  = internals[BETA_MIN];
        const auto& beta_max_constr  = internals[BETA_MAX];
        const auto& p_out_min_constr = internals[P_OUT_MIN];
        const auto& p_out_max_constr = internals[P_OUT_MAX];

        bool pass_min = beta_min_constr.check(press_rate);
        bool pass_max = beta_max_constr.check(press_rate);

        if (pass_min & pass_max)
            return press_rate;

        double beta_min = internals[BETA_MIN].value();
        double beta_max = internals[BETA_MAX].value();

        double beta;
        if (!pass_min)
            beta = beta_min;
        else if (!pass_max)
            beta = beta_max;

        p_out = beta * p_in;

        if (!p_out_min_constr.check(p_out))
            p_out = p_out_min_constr.value();
        else if (!p_out_max_constr.check(p_out))
            p_out = p_out_max_constr.value();

        return p_out / p_in;
    }


    } //end namespace control
}//end namespace shimmer

// 0. Check activation warning
// 1. t is really time t? why -1? It is not it-1?
// 2. Use pipes as by_pass....not to heavy? Maybe better just add real stations?
// 3. Marco: In PWD, Ki = ZiTiR.. what about R? at nodes or at pipe? => use c2_nodes, c2_pipes? combination?
// 4. Marco: Maybe here I should recheck constraint to be sure new pwd < pwd_nominal;