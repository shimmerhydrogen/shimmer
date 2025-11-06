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


#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "infrastructure_graph.h"
#include "solver/incidence_matrix.h"
#include "solver/conservation_matrices.h"
#include "solver/quality_tracking_solver.h"
#include "solver/viscosity.h"
#include "verify_test.h"
#include "infra/launch_solver.h"

#include <filesystem>
namespace fs = std::filesystem;

using namespace shimmer;

std::vector<double>
make_ref_sol() 
{
    return std::vector<double>
    ({
       15914750,16700000,15679766.00702241,14601780.11385979,
       15432913.968529,15419809.49878528,15324955.34871516,
       15089515.4449182,14841467.36010754,15867548.4218792,
       15820133.08801724,15772879.01614502,15725938.73121835,
       16447229.84289386,16034243.13131681,15433336.72994541,
       15105127.16443368,16508268.66367633,16299462.60586477,
       16197394.45197337,16002627.67097471,15632673.36777537,
       15577254.04352749,15525585.56911791,15477530.44706138,
       184.7856628896391,182.2353186013897,179.5244107405513,
       175.5611853156477,171.0113242681625,114.128867828297,
       113.9422751753469,113.5746337434375,113.0255796174109,
       112.294192306184,320.5345159441952,322.0892792973021,
       319.8770374054209,317.8566130961824,314.6286249068031,
       160.0903234920585,160.7441355878735,161.5656351176511,
       161.7830210400305,161.73422388476,85.932977453961,
       82.67892720787167,79.54051653982562,76.49096557194422,
       73.49312305054791,-114.128867828297,-480.6248394362543,
       0,480,72,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    });
}


std::vector<double>
make_ref_molfrac()
{   
    return std::vector<double>({
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0.9509163434030115,0,0.04908365659698849,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.585655180141763e-17,
    0.9991241587297426,0,0.0008758412702574483,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.735460231865488e-17,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0.9972738066618451,0,0.002726193338154802,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.07399264676449e-17,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0.9999337741924093,0,6.622580759089483e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.634176980975789e-17,
    1.0027791697896,0,-0.002779169789599525,0,0,0,0,0,0,0,0,0,0,-9.679584135589893e-20,0,0,0,0,0,0,-2.961075161702951e-17,
    0.4280621581358214,0,0.5719378418641785,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.829915861210626e-17,
    0.8119564815562101,0,0.1880435184437898,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.805252389741335e-17,
    0.9667721703182779,0,0.03322782968172203,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.298919487550558e-17,
    0.9972773792834427,0,0.002722620716557344,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.039351446021997e-18,
    0.02983427185903828,0,0.9701657281409617,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.559729905738403e-16,
    0.2762503719759692,0,0.7237496280240308,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.266341093346125e-17,
    0.6827632058777569,0,0.3172367941222431,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0.860376765166475,0,0.1396232348335251,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.738720904417735e-17,
    0.9992348097766004,0,0.0007651902233995892,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.026713639206805e-17,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    });
}

int main()
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
 
    shimmer::config cfg;
    cfg.refine = true;
    cfg.dx = 0;
    cfg.steps = 10;
    cfg.dt = 8100;
    cfg.tol_std = 1.e-4;
    cfg.quality_tracking = true;
    cfg.qt_steady = false;
    cfg.database = "../unit_tests/testdata/test_qt.db";
  
    //---------------------------------------------------------------    
    shimmer::infrastructure infra(cfg);
    shimmer::prepare_infrastructure(cfg, infra);
    //---------------------------------------------------------------
    shimmer::variable guess = initial_guess(infra);
    
    // Solver
    using solver_t = shimmer::qt_solver<shimmer::gerg_aga,
        shimmer::viscosity_type::Constant>;

    solver_t qt(infra, cfg.temperature, cfg.refine, cfg.steps);

    if(cfg.qt_steady)
        qt.steady_state(guess, cfg.dt_std, cfg.tol_std);  
    else
    {
        qt.initialization(guess, cfg.dt_std, cfg.tol_std);  
        qt.advance(cfg.dt, cfg.tol);
    }

    using vector_t =  Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using matrix_t =  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    vector_t sol_unstd  = qt.solution();
    matrix_t mlf  = qt.molar_fractions_evolution()[cfg.steps-1].transpose();
    vector_t mf_unstd = Eigen::Map<vector_t>(mlf.data(),
                                             mlf.size()); 

    auto sol_ref = make_ref_sol();
    auto mf_ref  = make_ref_molfrac();

    bool pass_sol =  verify_test("QT solver: solution ", sol_unstd, sol_ref); 
    bool pass_mf  =  verify_test("QT solver: molarfrac", mf_unstd, mf_ref); 

    bool pass = pass_sol && pass_mf;
    std::cout << pass  << std::endl; 

    return !(pass);

}
