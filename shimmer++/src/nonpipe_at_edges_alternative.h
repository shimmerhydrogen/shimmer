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
    };


    enum variable_to_control
    {
        PRESSURE_OUT,
        PRESSURE_IN,
        PIN_LOWER_POUT,
        BETA,
        FLUX,
        NONE,
    };


    enum control_type
    {
        BY_PASS,
        SHUT_OFF,
        POWER_DRIVER,
        PRESSURE_OUT,
        PRESSURE_IN,
        BETA,
        FLUX,
    };


    class constraint_control
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

        public:
        model(control_type ctype)
        {
            switch (ctype)
            {
            case control_type::BY_PASS:
                coeffs[0] =-1.0;
                coeffs[1] = 1.0;
                coeffs[2] = 0.0;
                coeffs[3] = 0.0;
                break;
            case control_type::SHUT_OFF:
                coeffs[0] = 0.0;
                coeffs[1] = 0.0;
                coeffs[2] = 1.0;
                coeffs[3] = 0.0;
                break;
            case control_type::POWER_DRIVER:
                free_index = std : vector<int>{ 0,1,2,3 };
                break;
            case control_type::PRESSURE_IN:
                coeffs[0] = 1.0;
                coeffs[1] = 0.0;
                coeffs[2] = 0.0;
                free_index = std:vector<int>{ 3 };
                break;
            case control_type::PRESSURE_OUT:
                coeffs[0] = 0.0;
                coeffs[1] = 1.0;
                coeffs[2] = 0.0;
                free_index = std:vector<int>{ 3 };
                break;
            case control_type::FLUX:
                coeffs[0] = 0.0;
                coeffs[1] = 0.0;
                coeffs[2] = 1.0;
                free_index = std:vector<int>{ 3 };
                break;
            default:
                break;
            }
        }

        set_coefficient(size_t index, double value)
        {
            assert(free_index.size() > 0);
            assert("Error: In model of edge station, coefficient at index is not modifiable" && std::find(free_index.begin(), free_index.end(), item) != vec.end())

            coeffs[index] = value;
        }
    };


    class control_state
    {
    public:
        model m_;
        constraint_control internal_;
        control_type type_;

        control_state(const model& m, const constraint_control& internal):
            m_(m), internal_(internal)
            {};
        control_state(const control_type& type, const constraint_control& internal):
            m_(model(type)), internal_(internal)
            {};

        virtual control_hard(double time = 0);

        set_c1(double value){ m.set_cofficient(0, value)};
        set_c2(double value){ m.set_cofficient(1, value)};
        set_c3(double value){ m.set_cofficient(2, value)};
        set_rhs(double value)
        {
            m.set_cofficient(3, value)
        };


    };


/*
    class shut_off : control_state
    { public:
        shut_off()
        {
            type_ = control_type::SHUT_OFF;
        };
        bool control_hard(double & p_in,
                          double & p_out,
                          double & G,
                          double & L,
                          const std::vector<>& internals)
        {
            return true; //p_out > p_in; //hard_.check(p_in, p_out);
        }
    };


    class by_pass: control_state
    {
        public:
        by_pass()
        {
            type_ = control_type::BY_PASS;
        };

        bool control_hard(double value, const std::vector<>& internals)
        {
            //p_in > p_out;
            return true;  //hard_.check(vel);
        }
    };
*/

    class control_power_driver: public control_state
    {
        double ramp_coeff_;

        public:

        power_driver(double ramp, const constraint& hard):
            type_(control_type::POWER_DRIVER),
            control_state(type_, hard),
            ramp_coeff_(ramp)
        {};


        bool control_hard(double time)
        {
            // 1. Check constraint
            bool pass = internal_.check(model.d, time);

            // 2. Control variable
            auto pw_nominal = internal_.value();

            if (pass)
                model.d = model.d * (t - 1.0) + ramp_coeff_ * pw_nominal;
                // WK: Maybe here I should recheck constraint to be sure new pwd < pwd_nominal;
            else
                model.d = pw_nominal;

            return pass;
        }
    };


    class control_rhs : public control_state
    {
        public:

        control_rhs(const control_type& type,
                    const constraint   & internal):
                    control_state(type, internal),
                    type_(type),
                    variable_name_(var_name)
        {
            assert("ERROR: Your control state allows modifications of fixed coefficients"
                && (type_ == control_type::BY_PASS ||
                    type_ == control_type::SHUT_OFF ||
                    type_ == control_type::FLUX ||
                    type_ == control_type::PRESSURE_IN ||
                    type_ == control_type::PRESSURE_OUT));
        };

        bool control_hard(double time)
        {
            bool pass = internal_.check(model.d);

            // Check constraint
            if (pass)
                return true;

            // Control variable
            model.d = internal_.value();

            return false;
        }
    };

/*
    class pressure_out : control_state
    {
            public:
        pressure_out()
        {
            type_ = control_type::PRESSURE_OUT;
        };

        bool control_hard(double & p_in,
                          double & p_out,
                          double & G,
                          double & L,
                          const std::vector<>& internals)
        {
            bool pass = hard_.check(p_out);

            // Check constraint
            if (pass)
            {
                model.d = p_out;
                return true;
            }

            // Control variable
            model.d = hard_.value();


            return false;
        }
    };


    class pressure_in: control_state
    {
        public:
        pressure_in()
        {
            type_ = control_type::PRESSURE_IN;
        };

        bool control_hard(double & p_in,
                          double & p_out,
                          double & G,
                          double & L,
                          const std::vector<>& internals)
        {
            bool pass = hard_.check(p_in);

            // Check constraint
            if (pass)
            {
                model.d = p_in;
                return true;
            }

            // Control variable
            model.d = hard_.value();

        }
    };


    class flux: control_state
    {
        public:
        flux()
        {
            type_ = control_type::FLUX;
        };


        bool control_hard(double & p_in,
                          double & p_out,
                          double & G,
                          double & L,
                          const std::vector<>& internals)
        {
            bool pass = hard_.check(G);

            // Check constraint
            if (pass)
            {
                model.d = G;
                return true;
            }

            // Control variable
            model.d = hard_.value();

            return false;
        }
    };

*/
    auto
    make_power_driver_control(double PWD_nominal, double ramp)
    {
        auto internal = constraint_control(hardness_type::HARD_CONTROL,
                                           constraint_type::GREATER_EQUAL,
                                           PWD_nominal);

        return control_power_driver(ramp, internal);
    }

    auto
    make_pressure_out_control(double pressure_out_max)
    {
        auto internal = constraint_control(hardness_type::HARD_CONTROL,
                    constraint_type::LOWER_EQUAL, pressure_out_max);

        return control_rhs(control_type::PRESSURE_OUT, internal, variable_to_control::PRESSURE_OUT);
    }

    auto
    make_pressure_in_control(double pressure_in_min)
    {
        auto internal = constraint_control(hardness_type::HARD_CONTROL,
                                           constraint_type::GREATER_EQUAL,
                                           pressure_in_min);

        return control_rhs(control_type::PRESSURE_IN, internal, variable_to_control::PRESSURE_IN);
    }

    auto
    make_by_pass_control(const variable_to_control & control_var)
    {
        auto internal = constraint_edge(hardness_type::HARD_CONTROL,
                                        constraint_type::GENERAL_EQUAL, 0);

        auto c = control_rhs(control_type::BY_PASS, internal, control_var);

        c.set_rhs(0);

        return c;
    }

    auto
    make_shutoff_control(const variable_to_control& control_var)
    {
        auto internal = constraint_edge(hardness_type::HARD_CONTROL,
                                        constraint_type::GENERAL_EQUAL, 1);

        auto c = control_rhs(control_type::SHUT_OFF, internal, control_var);

        c.set_rhs(0);

        return c;
    }
} //end namespace control

    class station
    {
        size_t count;
        std::string name_;
        std::vector<bool> active_;

        std::vector<constraint_edge> internals_;
        std::vector<constraint_edge> externals_;

        std::vector<control_state> controls_off;
        std::vector<control_state> controls_on;

        station(const std::string& name, const std::vector<bool>& active,
                     const std::vector<constraint_edge>& internals,
                     const std::vector<constraint_edge>& externals):
                     name_(name), active_(active),
                     internals_(internals), externals_(externals)
        {
            count = 0;
        };

        inline active(int step){ return active_[step]};

        inline which_control_type()
        {
            size_t idx = count % controls.size();
            controls[idx].type;
            return;
        };

        inline num_controls(){return controls.size();};

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


    template<typename CONSTR_NONPIPE>
    std::vector<CONSTR_NONPIPE>
    build_multiple_constraints(const std::vector<pair_input_t>&limits, hardness_type ht)
    {
       std::vector<CONSTR_NONPIPE> constr_vector(limits.size());
        size_t i = 0;
        for(const auto&  ec : user_limits)
            constr_vector[i++] = constraint(ht,
                                            limits[i].first,
                                            limits[i].second);
        return constr_vector;
    }


    auto
    make_valve(const std::vector<bool>& activate_history,
               double velocity_limit,
               const std::vector<pair_input_t>& user_limits)
    {
        auto internal = constraint_edge(hardness_type::HARD,
                                        constraint_type::GREATER_EQUAL,
                                         velocity_limit);
        auto externals = build_multiple_constraints<constraint_edge>(user_limits);

        station valve("VALVE", activate_history, internal, externals);

        auto c_by_pass = make_by_pass_control(variable_to_control::NONE);
        auto c_shutoff = make_shutoff_control(variable_to_control::NONE);

        valve.set_control(c_by_pass);
        valve.set_control(c_shutoff);

        return;
    }


    auto
    make_compressor(size_t control_node,
                    double ramp,
                    double efficiency,
                    double power_nominal,
                    const std::vector<double>& activate_history,
                    const std::vector<pair_input_t>& pressure_limits,
                    const std::vector<pair_input_t>& hard_limits,
                    const std::vector<pair_input_t>& user_limits)
    {
        auto thresholds = build_multiple_constraints<constraint>(pressure_limits, hardness_type::HARD);
        auto internals  = build_multiple_constraints<constraint_edge>(hard_limits, hardness_type::HARD);
        auto externals  = build_multiple_constraints<constraint_edge>(user_limits, hardness_type::SOFT);

        station compressor("COMPRESSOR", activate_history, internals, externals);

        compressor.set_control_node(control_node, thresholds);
        compressor.set_control_modes(std::vector<control_type>{POWER_DRIVER, P_OUT, P_IN, BETA, FLUX});

        auto c_by_pass = make_by_pass_control(variable_to_control::PIN_LOWER_POUT);
        auto c_shutoff = make_shutoff_control(variable_to_control::PIN_LOWER_POUT);

        compressor.set_control_off(c_by_pass);
        compressor.set_control_off(c_shutoff);

        return compressor;
    }


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

} //end namespace control
}//end namespace shimmer

// 1. t is really time t? why -1? It is not it-1?
// 2. Use pipes as by_pass....not to heavy? Maybe better just add real stations?
// 3. Marco: In PWD, Ki = ZiTiR.. what about R? at nodes or at pipe? => use c2_nodes, c2_pipes? combination?
// 4. Marco: Maybe here I should recheck constraint to be sure new pwd < pwd_nominal;