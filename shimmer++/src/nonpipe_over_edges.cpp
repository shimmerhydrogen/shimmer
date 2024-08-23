#include "../src/nonpipe_over_edges.h"

namespace shimmer
{
namespace edge_station
{

bool
control::constraint::check(double variable) const
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
        free_index = std::vector<int>{ 0,1,2,3 };
        control_index_ = 3;
        break;
    case control::type::PRESSURE_IN:
        coeffs[0] = 1.0;
        coeffs[1] = 0.0;
        coeffs[2] = 0.0;
        free_index = std::vector<int>{ 3 };
        control_index_ = 3;
        break;
    case control::type::PRESSURE_OUT:
        coeffs[0] = 0.0;
        coeffs[1] = 1.0;
        coeffs[2] = 0.0;
        free_index = std::vector<int>{ 3 };
        control_index_ = 3;
        break;
    case control::type::FLUX:
        coeffs[0] = 0.0;
        coeffs[1] = 0.0;
        coeffs[2] = 1.0;
        free_index = std::vector<int>{ 3 };
        control_index_ = 3;
        break;
    case control::type::BETA:
        coeffs[1] = 1.0;
        coeffs[2] = 0.0;
        coeffs[3] = 0.0;
        free_index = std::vector<int>{0};
        control_index_ = 1;
        break;
    default:
        break;
    }
}

void
control::model::set_coefficient(size_t index, double value)
{
    assert(free_index.size() > 0);
    assert("Error: In model of edge station, coefficient at index is not modifiable"
            && std::find(free_index.begin(), free_index.end(), index) != free_index.end());

    coeffs[index] = value;
}

void
control::model::set_control_coefficient(double value)
{
    coeffs[control_index_] = value;
}

double
control::model::coefficient(size_t index)
{
    assert(index >= 0 || index <= 3);

    return coeffs[index];
}

double
control::model::control_coefficient()
{
    return coefficient(control_index_);
}

bool
control::state::control_hard(size_t step_)
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
    auto internal = control::constraint(control::hardness_type::HARD,
                                        control::constraint_type::GREATER_EQUAL,
                                        PWD_nominal);
    return state(control::type::POWER_DRIVER, internal);
    //return control_power_driver(ramp, internal);
}

auto
control::make_pressure_out_control(double pressure_out_max)
{
    auto internal = control::constraint(hardness_type::HARD,
                constraint_type::LOWER_EQUAL, pressure_out_max);

    return state(control::type::PRESSURE_OUT, internal);
}

auto
control::make_pressure_in_control(double pressure_in_min)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                        constraint_type::GREATER_EQUAL,
                                        pressure_in_min);

    return state(control::type::PRESSURE_IN, internal);
}

auto
control::make_by_pass_control(const constraint_type& ctype)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    ctype, 0);

    return state(control::type::BY_PASS, internal);
}

auto
control::make_shutoff_control(const constraint_type& ctype)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    ctype, 1);

    return  state(control::type::SHUT_OFF, internal);
}

auto
control::make_beta_min_control(double beta_min)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    constraint_type::LOWER_EQUAL, beta_min);

    return  state(control::type::BETA, internal);
}

auto
control::make_beta_max_control(double beta_max)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    constraint_type::GREATER_EQUAL, beta_max);

    return  state(control::type::BETA, internal);
}

auto
control::make_flux_control(double flux_max)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    constraint_type::GREATER_EQUAL, flux_max);

    return  state(control::type::FLUX, internal);
}
//==========================================================

bool
station::control_hard(size_t step)
{
    //size_t idx = count % controls[idx].size();
    bool pass = mode_.control_hard(step);

    if (!pass)
    {
        count++;
        return false;
    }

    return true;
}

/*
bool
station::check_soft(double p, double l, size_t step)
{
    bool success = true;
    for(auto& e : externals_)
    {
        if(!e.check(p, l, step))
        {
            success = false;
            std::cout << "WARNING SOFT:" << name_  << " constraint violated." << std::endl;
            std::cout << " * Soft constraint (-) : " << e.type() << " " << e.value(step) <<std::endl;
            std::cout << " * (press , lrate) :  ("<<p << "," << l<< ") " << std::endl;
        }
    }

    return success;
}
*/

void
station::set_c1(double value)
{
    //size_t idx = count % mode_.size();
    mode_.model_.set_coefficient(0, value);
}

void
station::set_c2(double value)
{
    //size_t idx = count % mode_.size();
    mode_.model_.set_coefficient(1, value);
}

void
station::set_c3(double value)
{
    //size_t idx = count % controls.size();
    mode_.model_.set_coefficient(2, value);
}

void
station::set_rhs(double value)
{
    //size_t idx = count % controls.size();
    mode_.model_.set_coefficient(3, value);
};



void
station::activate(size_t step, double target_pressure)
{
    bool on   = active_history_[step];

    if (on)
        itor = controls_on.begin();
    else
        itor = controls_off.begin();

    return;
}


compressor::compressor(const std::string& name,
        double efficiency,
        double ramp_coeff,
        const std::vector<bool>& activate_history,
        const std::vector<control::constraint>& internals,
        const std::vector<control::constraint>& externals):
        station(name, activate_history, internals, externals)
{
    /*
    std::vector<external_type> mandatory_exts = { BETA_MIN,
    BETA_MAX,
    P_OUT_MAX,
    P_OUT_MIN,
    P_IN_MIN,
    P_THRESHOLD_MIN,
    P_THRESHOLD_MAX};

    for(const auto& e : mandatory_exts)
    {
        if (externals.find(e) == externals.end())
        {
            auto messsage = "Mandatory limit was not goven for" + e.key();
            throw std::exception(messsage);
        }
    }
    */
};


bool
compressor::control_hard(size_t step)
{

    //size_t idx = count % controls.size();
    bool pass = mode_.control_hard(step);

    if (!pass)
    {
        count++;
        return false;
    }

    if(mode_.type_ == control::type::POWER_DRIVER)
    {
        auto pwd = mode_.model_.control_coefficient();
        auto pwd_nominal = mode_.internal_.value();
        auto pwd_control = pwd * (step - 1.0) + ramp_coeff_ * pwd_nominal;
        mode_.model_.set_control_coefficient(pwd_control);
    }

    return true;
}

void
compressor::activate(size_t step, double target_pressure)
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

    size_t p_threshold_min_ind = 0;
    size_t p_threshold_max_ind = 1;

    double pn = target_pressure;
    bool pass_down = externals_[p_threshold_min_ind].check(pn);
    bool pass_up   = externals_[p_threshold_max_ind].check(pn);

    bool on   = (!active_history_[step] & pass_up)
                || (active_history_[step] & !pass_down);

    if (on)
    {
        itor_ = controls_on.begin();
        return;
    }

    if (p_out > p_in)
        itor_ = controls_off.begin() + 1;
    else
        itor_ = controls_off.begin() ;

    return;
}

double
compressor::compute_beta(double p_in,
                double p_out)
{
    size_t p_out_min_ind = 2;
    size_t p_out_max_ind = 2;
    size_t beta_min_ind = 3;
    size_t beta_max_ind = 4;

    assert(controls_on[p_out_max_ind].type_== external_type::P_OUT_MAX
        || controls_on[beta_max_ind].type_ == external_type::BETA_MAX
        || controls_on[beta_min_ind].type_ == external_type::BETA_MIN
        || externals_[p_out_min_ind].type() == external_type::P_OUT_MIN
    );

    double press_rate = p_out / p_in;

    const auto& beta_min_constr  = controls_on[beta_min_ind].internal_;
    const auto& beta_max_constr  = controls_on[beta_max_ind].internal_;
    const auto& p_out_max_constr = controls_on[p_out_max_ind].internal_;
    const auto& p_out_min_constr = controls_on[p_out_min_ind].internal_;

    bool pass_min = beta_min_constr.check(press_rate);
    bool pass_max = beta_max_constr.check(press_rate);

    if (pass_min & pass_max)
        return press_rate;

    double beta_min = beta_min_constr.value();
    double beta_max = beta_max_constr.value();

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


std::vector<control::constraint>
build_multiple_constraints(const std::vector<std::pair<control::constraint_type, double>>& limits,
                            control::hardness_type ht)
{
    std::vector<control::constraint> constr_vector(limits.size());
    size_t i = 0;
    for(const auto&  l : limits)
        constr_vector[i++] = control::constraint(ht,
                                        limits[i].first,
                                        limits[i].second);
    return constr_vector;
}


auto
make_regulator(double velocity_limit,
                const std::vector<bool>& activate_history,
                const std::vector<  std::tuple<external_type,
                                    control::constraint_type,
                                    double>> & user_limits)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                        control::constraint_type::GREATER_EQUAL,
                                        velocity_limit);
    auto externals = build_multiple_constraints(user_limits, control::hardness_type::SOFT);

    station regulator("REGULATOR", activate_history, internal, externals);

    auto c_by_pass = make_by_pass_control(control::constraint_type::GREATER_EQUAL);
    auto c_shutoff = make_shutoff_control(control::constraint_type::LOWER);

    regulator.add_control_on(c_by_pass);
    regulator.add_control_off(c_shutoff);

    return;
}

auto
make_valve(const std::vector<bool>& activate_history,
            double velocity_limit,
            const std::vector<  std::tuple<external_type,
                                control::constraint_type,
                                double>> & user_limits)
    {
    auto internal = control::constraint(control::hardness_type::HARD,
                                        control::constraint_type::GREATER_EQUAL,
                                        velocity_limit);
    auto externals = build_multiple_constraints(user_limits, control::hardness_type::SOFT);

    station valve("VALVE", activate_history, internal, externals);

    auto c_by_pass = make_by_pass_control(control::constraint_type::NONE);
    auto c_shutoff = make_shutoff_control(control::constraint_type::NONE);

    valve.add_control_on(c_by_pass);
    valve.add_control_off(c_shutoff);

    return;
}


auto
make_compressor(double ramp,
                double efficiency,
                const std::vector<double>& activate_history,
                const std::unordered_map<external_type,
                                        std::pair<control::constraint_type,
                                        double>> & user_limits)
{
    auto flux_hard = control::constraint(control::hardness_type::HARD,
                                         control::constraint_type::GREATER_EQUAL,
                                         0.0);

    std::vector<control::constraint> internals = { flux_hard };

    auto externals = build_multiple_constraints(
                                      { user_limits[P_THRESHOLD_MIN],
                                        user_limits[P_THRESHOLD_MAX],
                                        user_limits[P_OUT_MIN]
                                       }, control::hardness_type::SOFT);

    compressor cmp("COMPRESSOR", activate_history, internals, externals);

    auto c_by_pass = control::make_by_pass_control(control::constraint_type::EQUAL);
    auto c_shutoff = control::make_shutoff_control(control::constraint_type::EQUAL);

    cmp.add_control_off(c_by_pass);
    cmp.add_control_off(c_shutoff);

    auto c_power_driver = control::make_power_driver_control(user_limits[PWD_NOMINAL].second, ramp);
    auto c_press_in  = control::make_pressure_in_control(user_limits[P_IN_MIN].second);
    auto c_press_out = control::make_pressure_out_control(user_limits[P_OUT_MAX].second);
    auto c_beta_min  = control::make_beta_min_control(user_limits[BETA_MIN].second);
    auto c_beta_max  = control::make_beta_max_control(user_limits[BETA_MAX].second);
    auto c_flux = control::make_flux_control(user_limits[FLUX_MAX].second);

    cmp.add_control_on(c_power_driver);
    cmp.add_control_on(c_press_in);
    cmp.add_control_on(c_press_out);
    cmp.add_control_on(c_beta_min);
    cmp.add_control_on(c_beta_max);
    cmp.add_control_on(c_flux);

    return cmp;
}




} //end namespace control
}//end namespace shimmer
