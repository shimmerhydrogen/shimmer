#pragma once
#include <Eigen/Dense>
#include <memory>
#include <iostream>
#include <vector>
#include <cassert>

namespace shimmer
{

    enum control_type
    {
        BY_PASS,
        SHUT_OFF,
        POWER_DRIVER,
        P_OUT,
        P_IN,
        BETA,
        FLUX,
    };

    enum constraint_pipe_type
    {
        B_IN_MIN_GREATER_EQUAL,
        B_MAX_LOWER_EQUAL,
        P_OUT_MAX_LOWER_EQUAL,
        P_OUT_MIN_GREATER_EQUAL,
        P_IN_MIN_GREATER_EQUAL,
        P_IN_MAX_LOWER_EQUAL,
        FLUX_MIN_GREATER_EQUAL,
        FLUX_MAX_LOWER_EQUAL,
        V_MAX_LOWER_EQUAL,
        V_MIN_GREATER_EQUAL,
        PW_NOMINAL_LOWER_EQUAL,
    }


    class model
    {
        double c1;
        double c2;
        double c3;
        double d;
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
    make_valve(const std::vector<pair_input_t>& hard_limits,
               const std::vector<pair_input_t>& user_limits)
    {
        station_pipe valve("VALVE", internals, externals);
        valve.set_control_modes(std::vector<control_type>{BY_PASS, SHUT_OFF});

        return;
    }

    auto
    make_regulator(const std::vector<pair_input_t>& hard_limits,
                   const std::vector<pair_input_t>& user_limits)
    {
        station_pipe regulator("REGULATOR", internals, externals);
        regulator.set_control_modes(std::vector<control_type>{???});

        return;
    }

    auto
    make_compressor(size_t control_node,const std::vector<double>& activation_times,
            const std::vector<pair_input_t>& pressure_limits,
            const std::vector<pair_input_t>& hard_limits,
            const std::vector<pair_input_t>& user_limits)
    {
        auto thresholds = build_multiple_constraints<constraint>(pressure_limits, hardness_type::HARD);
        auto internals  = build_multiple_constraints<constraint_pipe>(hard_limits, hardness_type::HARD);
        auto externals  = build_multiple_constraints<constraint_pipe>(user_limits, hardness_type::SOFT);

        station_pipe compressor("COMPRESSOR", internals, externals);
        compressor.set_control_node(control_node, thresholds);
        compressor.set_control_modes(std::vector<control_type>{POWER_DRIVER, P_OUT, P_IN, BETA, FLUX});

        compressor.activate(activation_times);
        return compressor;
    }



    void
    compressor_activation()
    {
        if (active)
        {
            if (Pn < P_THRESHOLD_DOWN) //threshold_pressure_[0].check(Pn, 0, 0);
                return control_modes_[i];
        }
        else
        {
            if (Pn > P_THRESHOLD_UP) //threshold_pressure_[1].check(Pn, 0, 0);
                return control_modes_[i];
        }

        if(p_out > p_in)
            return control_type::SHUT_OFF;
        else
            return control_type::BY_PASS;
    }




    double
    check_beta()
    {
        double beta;

        if (pout < beta_min * p_in)
        {
            beta = beta_min;
            p_out = beta_min * p_in;

            if ()
            {
                p_out = p_out_min;
                beta = p_out/p_in;
            }
            if ()
            {
                p_out = p_out_max;
                beta = p_out/p_in;
            }
        }
        else if (p_out < beta_max * p_in)
        {
            beta = beta_max;
            p_out = beta_max * p_in;

            if ()
            {
                p_out = p_out_min;
                beta = p_out/p_in;
            }
            if ()
            {
                p_out = p_out_max;
                beta = p_out/p_in;
            }

        }
        else
            beta = p_out/p_in;

        return beta;
    }


    model
    control_power_driver(double t, double ramp, double p_out, double mass_flow)
    {
        auto beta = check_beta();
        internals[POWER_DRIVER].check();

        if (pw_driver < pw_nominal)
            pw_driver = pw_driver * (t - 1) + ramp * (pw_nominal);
        else
            pw_driver = pw_nominal;

        model m;
        model.c1 = -Ki * mass_flow * beta / p_in;
        model.c2 =  Ki * mass_flow * beta / p_out;
        model.c3 = (Ki/ck) * (beta - 1.0);
        model.d  =  pw_driver;

        return m;
    }


    auto
    control_pressure_in(double p_in)
    {
        model m;
        model.c1 = 1.0;
        model.c2 = 0.0;
        model.c3 = 0.0;
        model.d = p_in;
        return m;
    }


    auto
    control_pressure_out(double p_out)
    {
        model m;
        model.c1 = 0.0;
        model.c2 = 1.0;
        model.c3 = 0.0;
        model.d = p_out;
        return m;
    }

    auto control_flow_rate(double flow_rate)
    {
        model m;
        model.c1 = 0.0;
        model.c2 = 0.0;
        model.c3 = 1.0;
        model.d = flow_rate;

        return m;
    }

    auto
    control_shut_off()
    {
        model m;

        model.c1 = 0.0;
        model.c2 = 0.0;
        model.c3 = 1.0;
        model.d = 0.0;

        return m;
    }

    auto
    control_by_pass()
    {
        model m;

        model.c1 =-1.0;
        model.c2 = 1.0;
        model.c3 = 0.0;
        model.d = 0.0;

        return m;
    }


    class edge_station
    {
        bool active_;
        std::string name_;
        model m_;

        size_t control_index_;
        size_t threshold_pressure_node_;
        std::vector<control_type> control_modes_;

        std::vector<constraint> threshold_pressure_;
        std::vector<constraint_pipe> internals_;
        std::vector<constraint_pipe> externals_;


        station_pipe(const std::string& name,
            const std::vector<constraint_pipe>& internals,
            const std::vector<constraint_pipe>& externals_):
            name_(name), internals_(internals), externals_(externals)
        {};

        void
        add_control_node(size_t node_index, const std::vector<constraint>& pressure_limits)
        {
            threshold_pressure_ = pressure_limits;
            threshold_pressure_node_ = node_index;

            // TODO: Check node exists!
            return;
        }

        model
        control(double pressure_in, double pressure_out, double flux)
        {
            activation();

            switch (control_modes_[control_index_])
            {
            case control_type::SHUT_OFF:
                m = control_shut_off();
                break;
            case control_type::BY_PASS:
                m = control_by_pass();
                break;
            case control_type::POWER_DRIVER:
                m = control_power_driver();
                break;
            case control_type::P_OUT:
                m = control_pressure_out(p_out);
                break;
            case control_type::P_IN:
                m = control_pressure_in(p_in);
                break;
            case control_type::BETA:
                m = control_pressure_ratio();
                break;
            case control_type::FLUX:
                m = control_flow_rate();
                break;
            default:
                break;
            }

            return m;
        }
    };

// 1. REGULATOR: control modes and actions?
// 2. COMPRESSOR: activation times as sense? Or is it needed only for regulator and valve?
//                how do i compute power driver? and coefficients Ki, ci?
// 3. NONPIPES_pipes iterate in the same cycle of NONPIPE_nodes?
} //end namespace shimmer