#include "../src/nonpipe_over_edges.h"

namespace shimmer
{
namespace edge_station
{

bool  control::constraint::check(double variable, size_t step)
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

control::model::model(control::type ctype)
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

void control::model::set_coefficient(size_t index, double value)
{
    assert(free_index.size() > 0);
    assert("Error: In model of edge station, coefficient at index is not modifiable" && std::find(free_index.begin(), free_index.end(), item) != vec.end())

    coeffs[index] = value;
}
void control::model::set_control_coefficient(double value)
{
    coeffs[control_index_] = value;
}
double control::model::control_coefficient()
{
    if (control_index_ < 0)
        return INF;

    return coeffs[control_index_];
}


virtual bool control::state::control_hard(double time)
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

auto
control::make_power_driver_control(double PWD_nominal, double ramp)
{
    auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                        constraint_type::GREATER_EQUAL,
                                        PWD_nominal);
    return state(control::type::POWER_DRIVER, internal);
    //return control_power_driver(ramp, internal);
}

auto
control::make_pressure_out_control(double pressure_out_max)
{
    auto internal = control::constraint(hardness_type::HARD_CONTROL,
                constraint_type::LOWER_EQUAL, pressure_out_max);

    return state(control::type::PRESSURE_OUT, internal);
}

auto
control::make_pressure_in_control(double pressure_in_min)
{
    auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                        constraint_type::GREATER_EQUAL,
                                        pressure_in_min);

    return state(control::type::PRESSURE_IN, internal);
}

auto
control::make_by_pass_control(const constraint_type& ctype)
{
    auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                    ctype, 0);

    return state(control::type::BY_PASS, internal);
}

auto
control::make_shutoff_control(const constraint_type& ctype)
{
    auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                    ctype, 1);

    return  state(control::type::SHUT_OFF, internal);
}

auto
control::make_beta_min_control(double beta_min)
{
    auto internal = control::constraint(hardness_type::HARD_CONTROL,
                                    constraint_type::LOWER_EQUAL, beta_min);

    return  state(control::type::BETA, internal);
}

auto
control::make_beta_max_control(double beta_max)
{
    auto internals = control::constraint(hardness_type::HARD_CONTROL,
                                    constraint_type::GREATER_EQUAL, beta_max);

    return  state(control::type::BETA, internal);
}

auto
control::make_flux_control(double flux_max)
{
    auto internals = control::constraint(hardness_type::HARD_CONTROL,
                                    constraint_type::GREATER_EQUAL, flux_max);

    return  state(control::type::FLUX, internal);
}
//==========================================================

bool
station::control_hard(double time)
{
    size_t idx = count % controls.size();
    bool pass = controls[idx].control_hard(time);

    if (!pass)
    {
        count++;
        return false;
    }
}

bool
station::check_soft(double p, double l, size_t step)
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

void
station::set_c1(double value)
{
    size_t idx = count % controls.size();
    controls[idx].set_coefficient(0, value);
}

void
station::set_c2(double value)
{
    size_t idx = count % controls.size();
    controls[idx].set_coefficient(1, value);
}

void
station::set_c3(double value)
{
    size_t idx = count % controls.size();
    controls[idx].set_coefficient(2, value);
}

void
station::set_rhs(double value)
{
    size_t idx = count % controls.size();
    controls[idx].set_coefficient(3, value);
};

bool
station::is_active(size_t step, double target_pressure)
{
    /* Warning: control node is not set for the moment.
    It is provided instead the pressure at the target node.
    The graph is not passed to the station, since the graph
    need station to be defined. Maybe improve by saving the
    control node and pass var.

    auto t = target(*itor, g);
    auto control_node = g[s].node_num;
    auto pn = var_.pressure[control_node];
    */

    double pn = target_pressure;
    bool pass_down = externals_[BETA_MIN].check(pn);
    bool pass_up   = externals_[BETA_MAX].check(pn);

    bool is_on   = (!active_history_[step] & pass_up)
                        || (active_history_[step] & !pass_down);

    return is_on;
}

bool
compressor::control_hard(size_t step, const infrastructure_graph& graph, const variable& var)
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


double
compressor_beta(double p_in,
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
