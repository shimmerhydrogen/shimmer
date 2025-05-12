/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../src/nonpipe_over_edges.h"

namespace shimmer
{
namespace edge_station
{

bool
control::constraint::check(double v) const
{
    switch (type_)
    {
        case LOWER_EQUAL:
            return (v <= value_);
        case GREATER_EQUAL:
            return (v >= value_);
        case EQUAL:
            return (std::abs((v - value_)/value_) < 1.e-14);
        case LOWER:
            return (v < value_);
        case GREATER:
            return (v > value_);
        case NONE:
            return true;
        default:
            throw std::invalid_argument("Boundary conditions not specified");
    }
};

control::model::model(compressor_mode ctype)
{
    coeffs = std::vector<double>{1.e20,1.e20,1.e20,1.e20};

    switch (ctype)
    {
    case compressor_mode::OFF_BYPASS:
        coeffs[0] =-1.0;
        coeffs[1] = 1.0;
        coeffs[2] = 0.0;
        coeffs[3] = 0.0;
        control_index_ = -1;
        break;
    case compressor_mode::OFF_CLOSED:
        coeffs[0] = 0.0;
        coeffs[1] = 0.0;
        coeffs[2] = 1.0;
        coeffs[3] = 0.0;
        control_index_ = -1;
        break;
    case compressor_mode::ON_POWER:
        free_index = std::vector<int>{ 0,1,2,3 };
        control_index_ = 3;
        break;
    case compressor_mode::ON_IPRESS:
        coeffs[0] = 1.0;
        coeffs[1] = 0.0;
        coeffs[2] = 0.0;
        free_index = std::vector<int>{ 3 };
        control_index_ = 3;
        break;
    case compressor_mode::ON_OPRESS:
        coeffs[0] = 0.0;
        coeffs[1] = 1.0;
        coeffs[2] = 0.0;
        free_index = std::vector<int>{ 3 };
        control_index_ = 3;
        break;
    case compressor_mode::ON_MASSFLOW:
        coeffs[0] =  0.0;
        coeffs[1] =  0.0;
        coeffs[2] =  1.0;
        free_index = std::vector<int>{ 3 };
        control_index_ = 3;
        break;
    case compressor_mode::ON_RATIO:
        coeffs[1] = 1.0;
        coeffs[2] = 0.0;
        coeffs[3] = 0.0;
        free_index = std::vector<int>{0};
        control_index_ = 1;
        break;
    default:
        std::cout << "ERROR: station does not know this control type => "<< +ctype << "\n";
        throw std::exception();
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
control::model::get_coefficient(size_t index) const
{
    assert(index >= 0 || index <= 3);

    return coeffs[index];
}

double
control::model::get_control_coefficient() const
{
    return get_coefficient(control_index_);
}

bool
control::mode::check_hard() const
{
    // 1. Get stored value in model
    auto value = model_.get_control_coefficient();

    // 2. Check constraint
    bool pass = internal_.check(value);

    return pass;
}

bool
control::mode::control_hard()
{
    // 1. Get stored value in model
    auto value = model_.get_control_coefficient();

    // 2. Check constraint
    bool pass = internal_.check(value);

    // 3. Control value with hard limit
    if (!pass)
        model_.set_control_coefficient(internal_.value());

    return pass;
}

auto
control::make_power_driver_mode(double PWD_nominal, double ramp)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                        control::constraint_type::GREATER_EQUAL,
        PWD_nominal);

    return mode(compressor_mode::ON_POWER, internal);
}

auto
control::make_pressure_out_mode(double pressure_out_max)
{
    auto internal = control::constraint(hardness_type::HARD,
                constraint_type::LOWER_EQUAL, pressure_out_max);

    return mode(compressor_mode::ON_OPRESS, internal);
}

auto
control::make_pressure_in_mode(double pressure_in_min)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                        constraint_type::GREATER_EQUAL,
                                        pressure_in_min);

    return mode(compressor_mode::ON_IPRESS, internal);
}

auto
control::make_bypass_mode(const constraint_type& ctype)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    ctype, 0);

    return mode(compressor_mode::OFF_BYPASS, internal);
}

auto
control::make_shutoff_mode(const constraint_type& ctype)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    ctype, 1);

    return  mode(compressor_mode::OFF_CLOSED, internal);
}


auto
control::make_beta_min_mode(double beta_min)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    constraint_type::LOWER_EQUAL, beta_min);

    return  mode(compressor_mode::ON_RATIO, internal);
}

auto
control::make_beta_max_mode(double beta_max)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    constraint_type::GREATER_EQUAL, beta_max);

    return  mode(compressor_mode::ON_RATIO, internal);
}

auto
control::make_flux_mode(double flux_max)
{
    auto internal = control::constraint(control::hardness_type::HARD,
                                    constraint_type::LOWER_EQUAL, flux_max);

    return  mode(compressor_mode::ON_MASSFLOW, internal);
}
//-----------------------------------------------------------------------------

bool
station::control_hard()
{
    return mode_.control_hard();
}

/*
bool
station::check_soft(double p, double l, size_t step)
{
    bool success = true;


    for(auto& e : control_on)
    for(size_t i = 0; i < control.size(); i++)
    {
        size_t idx = count % controls.size();


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
    mode_.model_.set_coefficient(0, value);
}

void
station::set_c2(double value)
{
    mode_.model_.set_coefficient(1, value);
}

void
station::set_c3(double value)
{
    mode_.model_.set_coefficient(2, value);
}

void
station::set_rhs(double value)
{
    mode_.model_.set_coefficient(3, value);
};


void
station::activate( size_t step,
                    int source_num,
                    int target_num,
                    const variable& var)
{
    on_   = active_history_[step];

    if (on_)
        mode_ = controls_on[0];
    else
        mode_ = controls_off[0];

    return;
}


void
edge_station::station::fill_model(control::mode& m,
                    int pipe_num,
                    int source_num,
                    int target_num,
                    const variable& var,
                    const vector_t& c2_nodes)
{
    auto p_in = var.pressure[source_num];
    auto p_out = var.pressure[target_num];

    switch (m.type_)
    {
        case compressor_mode::OFF_CLOSED:
        case compressor_mode::OFF_BYPASS:
            m.set_c3(p_in < p_out);
            break;
        case compressor_mode::ON_IPRESS:
            m.set_rhs(p_in);
            break;
        case compressor_mode::ON_OPRESS:
            m.set_rhs(p_out);
            break;
        case compressor_mode::ON_MASSFLOW:
            m.set_rhs(var.flux[pipe_num]);
            break;
        default:
            std::cout << "ERROR: station does not know this control type.\n";
            throw std::exception();
    }
}


std::ostream& operator<<(std::ostream& ofs, const station& st)
{
    ofs << " name : " << st.name_ << "\n";
    ofs << " current mode_type  : " << +st.which_mode_type() << "\n";

    return ofs;
}


//-----------------------------------------------------------------------------

compressor::compressor(const std::string& name,
                        double efficiency,
                        double ramp_coeff,
                        const std::vector<bool>& activate_history,
                        const std::vector<control::constraint>& internals,
                        const std::vector<control::constraint>& externals):
                        station(name, activate_history, internals, externals),
                        efficiency_(efficiency), ramp_coeff_(ramp_coeff)
{

    #if 0
    std::vector<external_type> mandatory_exts = { BETA_MIN,
                                                BETA_MAX,
                                                P_OUT_MAX,
                                                P_OUT_MIN,
                                                P_IN_MIN,
                                                P_THRESHOLD_MIN,
                                                P_THRESHOLD_MAX};

    for (const auto& e : mandatory_exts)
    {
        auto has_ext = [&](const control::constraint& c)
                        {
                            return c.type() == e;
                        };

        auto it = std::find(externals.begin(), externals.end(), has_ext);

        if (it == externals.end())
            throw std::invalid_argument("Mandatory limit was not given");
    }

    #endif
};


bool
compressor::control_hard()
{
    bool pass = mode_.control_hard();

    if (!pass)
        return false;

    if(mode_.type_ == compressor_mode::ON_POWER)
    {
        auto pwd = mode_.model_.get_control_coefficient();
        auto pwd_nominal = mode_.internal_.value();

        // TODO: This needs to be finished by adding pwd_old_ = pwd(t^{n-1})
        // auto pwd_control = pwd_old_ + ramp_coeff_ * pwd_nominal;
        // mode_.model_.set_control_coefficient(pwd_control);
    }

    return true;
}

/*
void
compressor::activate(size_t step,
                     int source_num,
                     int target_num,
                     const shimmer::variable& var)
{

    // WARNING: control_node set as default as target node
    auto control_num = target_num;

    const auto& nodes_pressure = var.pressure;
    double p_control = nodes_pressure[control_num];
    double p_in  = nodes_pressure[source_num];
    double p_out = nodes_pressure[target_num];

    size_t p_threshold_min_ind = 0;
    size_t p_threshold_max_ind = 1;

    bool pass_down = externals_[p_threshold_min_ind].check(p_control);
    bool pass_up   = externals_[p_threshold_max_ind].check(p_control);

    on_   = (!active_history_[step] & pass_up)
                || (active_history_[step] & !pass_down);

    if (on_)
    {
        mode_ = controls_on[0];
        return;
    }

    if (p_out > p_in)
    {
        assert(controls_off[0].type_ == control::mode_type::BY_PASS);
        mode_ = controls_off[0];
    }
    else
    {
        assert(controls_off[1].type_ == control::mode_type::SHUT_OFF);
        mode_ = controls_off[1];
    }

    return;
}
*/

double
compressor::compute_beta(double p_in,
                         double p_out)
{
    size_t beta_min_ind = 4;
    size_t beta_max_ind = 5;

    assert(externals_[beta_max_ind].type() == external_type::BETA_MAX
        || externals_[beta_min_ind].type() == external_type::BETA_MIN);

    const auto& beta_min_constr = externals_[beta_min_ind];
    const auto& beta_max_constr = externals_[beta_max_ind];

    double beta = p_out / p_in;

    bool pass_beta_min = beta_min_constr.check(beta);
    bool pass_beta_max = beta_max_constr.check(beta);

    double beta_min = beta_min_constr.value();
    double beta_max = beta_max_constr.value();

    if (!pass_beta_min)
        beta = beta_min;
    else if (!pass_beta_max)
        beta = beta_max;
    else
        return beta;

    size_t pout_min_ind = 2;
    size_t pout_max_ind = 3;

    assert(externals_[pout_max_ind].type() == external_type::P_OUT_MAX
        || externals_[pout_min_ind].type() == external_type::P_OUT_MIN);

    const auto& pout_max_constr = externals_[pout_max_ind];
    const auto& pout_min_constr = externals_[pout_min_ind];

    p_out = beta * p_in;

    bool pass_pout_min = pout_min_constr.check(p_out);
    bool pass_pout_max = pout_max_constr.check(p_out);

    bool pout_min = pout_min_constr.value();
    bool pout_max = pout_max_constr.value();

    if (!pass_pout_min)
        p_out = pout_min;
    else if (!pass_pout_max)
        p_out = pout_max;

    if ((!pass_pout_min && !pass_beta_min) || (!pass_pout_max && !pass_beta_max))
        std::cout << "WARNING: This condition shouldn't happend.";

    return p_out / p_in;
}


void
compressor::fill_model( control::mode& m,
                        int pipe_num,
                        int source_num,
                        int target_num,
                        const variable& var,
                        const vector_t& c2_pipes)
{

    auto p_in  = var.pressure[source_num];
    auto p_out = var.pressure[target_num];


    switch (m.type_)
    {
        case compressor_mode::OFF_CLOSED:
        case compressor_mode::OFF_BYPASS:
            m.set_c3(p_in < p_out);
            std::cout<< "I am inside fill_model with BYPASS or SHUTOFF \n"; 
            break;
        case compressor_mode::ON_RATIO:
        {
            auto beta = -beta_; // p_out / p_in;
            m.set_c1(beta);
            std::cout<< "I am inside fill_model with BETA \n"; 
            break;
        }
        case compressor_mode::ON_POWER:
        {
            double gamma = 1.4; // Or read from GERG
            double ck = (gamma - 1.0) / gamma;
            double beta =  p_out / p_in; //For tessting purposes, here we use simplified beta //compute_beta(p_in, p_out);
            double betaPowck = std::pow(beta, ck);
            double ZTR = c2_pipes[pipe_num];
            double K = ZTR / efficiency_;
            double G = var.flux[pipe_num];
            double KGB = K * G * betaPowck;

            double c1 = -KGB / p_in;
            double c2 =  KGB / p_out;
            double c3 = (K / ck) * (betaPowck - 1.0);
            double computed_pwd = c3 * G;

            m.set_c1(-KGB / p_in);
            m.set_c2( KGB / p_out);
            m.set_c3( c3);
            m.set_rhs(pwd_);//For testing purposes, here we use simplified pwd // computed_pwd);

            std::cout<< "I am inside fill_model with POWER_DRIVER \n"; 
            break;
        }
        case compressor_mode::ON_IPRESS:
            //m.set_rhs(p_in);
            m.set_rhs(m.internal_.value());
            break;
        case compressor_mode::ON_OPRESS:
            //m.set_rhs(p_out);
            m.set_rhs(m.internal_.value());
            std::cout<< "I am inside fill_model with PRESSURE_OUT \n"; 
            break;
        case compressor_mode::ON_MASSFLOW:
            //m.set_rhs(var.flux[pipe_num]);
            m.set_rhs(flux_);
            break;
        default:
            std::cout << "ERROR: compressor does not know this control type.\n";
            throw std::exception();
    }
}


std::ostream& operator<<(std::ostream& ofs, const compressor& cmp)
{
    ofs << " name : " << cmp.name_ << "\n";
    ofs << " ramp_coeff : " << cmp.ramp_coeff_ << "\n";
    ofs << " efficiency : " << cmp.efficiency_ << "\n";
    ofs << " current mode_type : " << +cmp.which_mode_type() << "\n";

    return ofs;
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
    auto hard_limit = control::constraint(  control::hardness_type::HARD,
                                            control::constraint_type::GREATER_EQUAL,
                                            velocity_limit);

    auto internal = std::vector<control::constraint>{ hard_limit };
    auto externals = build_multiple_constraints(user_limits, control::hardness_type::SOFT);

    station regulator("REGULATOR", activate_history, internal, externals);

    auto c_by_pass = make_bypass_mode(control::constraint_type::GREATER_EQUAL);
    auto c_shutoff = make_shutoff_mode(control::constraint_type::LOWER);

    regulator.add_mode_on(c_by_pass);
    regulator.add_mode_off(c_shutoff);

    return;
}


auto
make_valve(const std::vector<bool>& activate_history,
            double velocity_limit,
            const std::vector<  std::tuple<external_type,
                                control::constraint_type,
                                double>> & user_limits)
{
    auto hard_limit = control::constraint(  control::hardness_type::HARD,
                                            control::constraint_type::GREATER_EQUAL,
                                            velocity_limit);

    auto internal = std::vector<control::constraint>{ hard_limit };
    auto externals = build_multiple_constraints(user_limits, control::hardness_type::SOFT);

    station valve("VALVE", activate_history, internal, externals);

    auto c_by_pass = make_bypass_mode(control::constraint_type::NONE);
    auto c_shutoff = make_shutoff_mode(control::constraint_type::NONE);

    valve.add_mode_on(c_by_pass);
    valve.add_mode_off(c_shutoff);

    return;
}


compressor
make_compressor(double ramp,
                double efficiency,
                const std::vector<bool>& activate_history,
                const std::vector<std::pair<compressor_mode,double>>& modes_type_vec,
                std::unordered_map<external_type,
                                        std::pair<control::constraint_type,
                                        double>> & user_limits)                                                   
{
    auto flux_limit = control::constraint(control::hardness_type::HARD,
                                         control::constraint_type::GREATER_EQUAL,
                                         0.0);

    auto internals = std::vector<control::constraint>({ flux_limit });
    auto externals = build_multiple_constraints({user_limits[P_THRESHOLD_MIN],
                                                 user_limits[P_THRESHOLD_MAX],
                                                 user_limits[P_OUT_MIN],
                                                 user_limits[P_OUT_MAX],
                                                 user_limits[BETA_MIN],
                                                 user_limits[BETA_MAX]},
                                                control::hardness_type::SOFT);

    compressor cmp("COMPRESSOR", efficiency, ramp, activate_history, internals, externals);

    for(const auto & m : modes_type_vec)
    {
        auto model_name = m.first;
        auto model_value = m.second; 
        switch (model_name)
        {
            case compressor_mode::OFF_CLOSED:
            {
                std::cout << "Adding mode OFF................ SHUT_OFF \n";
                auto c_shutoff = control::make_shutoff_mode(control::constraint_type::EQUAL);
                cmp.add_mode_off(c_shutoff);
                break;
            }    
            case compressor_mode::OFF_BYPASS:
            {
                std::cout << "Adding mode OFF................ BY_PASS \n";
                auto c_by_pass = control::make_bypass_mode(control::constraint_type::EQUAL);
                cmp.add_mode_off(c_by_pass);
                break;
            }
            case compressor_mode::ON_IPRESS:
            {
                std::cout << "Adding mode ON ................ PRESS_IN \n";
                if(model_value != user_limits[P_IN_MIN].second)
                    throw std::invalid_argument("Minimun pressure for control != to user limits pressure_in_min");

                auto c_press_in  = control::make_pressure_in_mode(user_limits[P_IN_MIN].second);
                cmp.add_mode_on(c_press_in);
                break;
            }
            case compressor_mode::ON_OPRESS:
            {
                std::cout << "Adding mode ON ................ PRESS_OUT \n";
                if(model_value != user_limits[P_OUT_MAX].second)
                    throw std::invalid_argument("Maximun pressure for control != to user limits pressure_out_max");
                auto c_press_out = control::make_pressure_out_mode(user_limits[P_OUT_MAX].second);
                cmp.add_mode_on(c_press_out);
                break;
            }
            case compressor_mode::ON_RATIO:
            {
                std::cout << "Adding mode ON ................ BETA \n";
                
                cmp.beta_ = model_value;                
                auto c_beta_min  = control::make_beta_min_mode(user_limits[BETA_MIN].second);
                cmp.add_mode_on(c_beta_min);
                // This should add both controls on beta. But When I do so, there are two controls that call fill model,
                // so the rhs is fill twice.

                break;
            }
            case compressor_mode::ON_MASSFLOW:
            {
                std::cout << "Adding mode ON ................ FLUX \n";
                cmp.flux_ = model_value;
                auto c_flux = control::make_flux_mode(user_limits[FLUX_MAX].second);
                cmp.add_mode_on(c_flux);
                break;
            }
            case compressor_mode::ON_POWER:
            {
                std::cout << "Adding mode ON ................ POWER_DRIVER \n";
                cmp.pwd_ = model_value;
                auto c_power_driver = control::make_power_driver_mode(user_limits[PWD_NOMINAL].second, ramp);
                cmp.add_mode_on(c_power_driver);
                break;
            }
            default:
                std::cout << "ERROR: station does not know this control type.\n";
                throw std::exception();
        }

        // Init current mode with the first mode in the list of modes ON (bypass and shutoff does not count)
        cmp.change_mode_on(0);

        return cmp;
    }

    // Check these ones! It is seems there are not well defined. For example there ir also a control mode for BETA
    /*
    auto c_beta_min  = control::make_beta_min_mode(user_limits[BETA_MIN].second);
    cmp.add_mode_on(c_beta_min);

    auto c_beta_max  = control::make_beta_max_mode(user_limits[BETA_MAX].second);
    cmp.add_mode_on(c_beta_max);
    */
    return cmp;
}



    } //end namespace control
}//end namespace shimmer
