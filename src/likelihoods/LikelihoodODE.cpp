#include "Utils.h"
#include "LikelihoodODE.h"
#include "ODESolverCVODE.h"
#include "ProbabilityDistributions.h"
#include "CSVParser.h"
#include "boost/filesystem.hpp"   
#include <iostream>

LikelihoodODE::LikelihoodODE(size_t sampling_threads, size_t evaluation_threads)
{
}

LikelihoodODE::~LikelihoodODE()
{
}

bool LikelihoodODE::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
    size_t num_dynamic_variables = 8;
    size_t num_inference_variables = 29;

    this->varset = varset;
    if (varset->GetNumVariables() != num_inference_variables) {
        LOGERROR("Incorrect number of parameters");
        return false;
    }

    solver = std::make_shared<ODESolverCVODE>();

	ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodODE::CalculateDerivative, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4);
	solver->SetDerivativeFunction(derivative);
	solver->Initialize(num_dynamic_variables, NULL);
	solver->SetTolerance(1e-7, 1e-3);

    // Allocate parameter vector to number of parameters that are used
    parameter_values.resize(num_inference_variables);

    boost::filesystem::path cwd = boost::filesystem::current_path() / "normalized_oscillations.csv";
    
    bcm3::CSVParser parser;
    parser.Parse("/Users/huubvdent/Documents/Internship/BCM_RUNS/7d_runs_v7_3/example/normalized_oscillations.csv", ",", false);
    size_t num_timepoints = parser.GetNumColumns();
    timepoints.resize(num_timepoints);
    
    for(size_t i = 0; i < num_timepoints; i++){
        float x = stof(parser.GetColumnName(i));
        timepoints(i) = x;
    }


    // Timepoints at which the ODE solver should return values of the dynamic variables
    // Real max_time = 21600.0;
    // size_t num_timepoints = 216;
    // timepoints.resize(num_timepoints);
    // for (size_t i = 0; i < num_timepoints; i++) {
    //     timepoints(i) = 0 + max_time * (i / (Real)(num_timepoints - 1));
    // }

	return true;
}

bool LikelihoodODE::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
    // Apply variable transformations if those were specified in the prior
    for (ptrdiff_t i = 0; i < values.size(); i++) {
        parameter_values(i) = varset->TransformVariable(i, values(i));
    }

    // The initial conditions can be dependent on the parameters
    OdeVectorReal initial_conditions(8);
    initial_conditions(0) = parameter_values(21);
    initial_conditions(1) = parameter_values(22);
    initial_conditions(2) = parameter_values(23);
    initial_conditions(3) = parameter_values(24);
    initial_conditions(4) = parameter_values(25);
    initial_conditions(5) = parameter_values(26);
    initial_conditions(6) = parameter_values(27);
    initial_conditions(7) = parameter_values(28);

    boost::filesystem::path cwd = boost::filesystem::current_path() / "normalized_oscillations.csv";
    
    bcm3::CSVParser parser;
    parser.Parse("/Users/huubvdent/Documents/Internship/BCM_RUNS/7d_runs_v7_3/example/normalized_oscillations.csv", ",", false);

    // "/Users/huubvdent/Documents/Internship/BCM_RUNS/7d_runs_v7_3/example/normalized_oscillations.csv"

    // Integrate the ODE system and calculate the likelihood based on the solution
    // In this example we compare the first dynamic variable to a cosine
    if (solver->Simulate(initial_conditions.data(), timepoints, simulated_trajectories)) {
        logp = 0.0;
        for (size_t i = 0; i < timepoints.size(); i++) {
            // Real cosvalue = 100.0 * cos(timepoints(i) / 1140) + 150.0;
            Real datavalue = parser.GetEntry(0, i);
            logp += bcm3::LogPdfTnu3(datavalue, parameter_values(18) * simulated_trajectories(0, i) + parameter_values(19), parameter_values(20));
        }
    } else {
        logp = -std::numeric_limits<Real>::infinity();
    }

	return true;
}

// bool LikelihoodODE::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
// {
//     // Calculate the right-hand side of the differential equation
//     // The current timepoint and values of the dynamic variables are provided in t and y
//     // The derivative should be stored in dydt

//     Real mek_kcat_erk = parameter_values[0];
//     Real mek_km_erk = parameter_values[1];
//     Real vmax_erk_dephos = parameter_values[2];
//     Real km_erk_dephos = parameter_values[3];
//     Real erk_vmax_mek = parameter_values[4];
//     Real erk_km_mek = parameter_values[5];
//     Real erk_kI = parameter_values[6];
//     Real vmax_mek_dephos = parameter_values[7];
//     Real km_mek_dephos = parameter_values[8];
    
//     Real erk_pp = y[0];
//     Real mek_pp = y[1];
//     Real erk = y[2];
//     Real mek = y[3];

//     Real& derk_pp = dydt[0];
//     Real& dmek_pp = dydt[1];
//     Real& derk = dydt[2];
//     Real& dmek = dydt[3];
    
//     derk_pp = mek_kcat_erk * mek_pp * erk / (mek_km_erk + erk) - vmax_erk_dephos * erk_pp / (km_erk_dephos + erk_pp);
//     dmek_pp = erk_vmax_mek * mek / (erk_km_mek * (1 + erk_pp / erk_kI) + mek) - vmax_mek_dephos * mek_pp / (km_mek_dephos + mek_pp);
//     derk = vmax_erk_dephos * erk_pp / (km_erk_dephos + erk_pp) - mek_kcat_erk * mek_pp * erk / (mek_km_erk + erk);
//     dmek = vmax_mek_dephos * mek_pp / (km_mek_dephos + mek_pp) - erk_vmax_mek * mek / (erk_km_mek * (1 + erk_pp / erk_kI) + mek);

//     return true;
// }

// bool LikelihoodODE::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
// {
//     // Calculate the right-hand side of the differential equation
//     // The current timepoint and values of the dynamic variables are provided in t and y
//     // The derivative should be stored in dydt

//     Real erk_dephos = parameter_values[0];
//     Real erk_kI_raf = parameter_values[1];
//     Real erk_km_raf = parameter_values[2];
//     Real erk_vmax_raf = parameter_values[3];
//     Real mek_dephos = parameter_values[4];
//     Real mek_kcat_erk = parameter_values[5];
//     Real mek_km_erk = parameter_values[6];
//     Real raf_dephos = parameter_values[7];
//     Real raf_kcat_mek = parameter_values[8];
//     Real raf_km_mek = parameter_values[9];

//     Real erk_pp = y[0];
//     Real mek_pp = y[1];
//     Real raf_pp = y[2];
//     Real erk = y[3];
//     Real mek = y[4];
//     Real raf = y[5];

//     Real& derk_pp = dydt[0];
//     Real& dmek_pp = dydt[1];
//     Real& draf_pp = dydt[2];
//     Real& derk = dydt[3];
//     Real& dmek = dydt[4];
//     Real& draf = dydt[5];
    
//     derk_pp = mek_kcat_erk * mek_pp * erk / (erk + mek_km_erk) - erk_pp * erk_dephos;
//     dmek_pp = - mek_pp * mek_dephos + raf_kcat_mek * raf_pp * mek /(mek + raf_km_mek);
//     draf_pp = raf_dephos * raf_pp - erk_vmax_raf * raf / (erk_km_raf * (1 + erk_pp / erk_kI_raf) + raf);
//     derk = mek_kcat_erk * mek_pp * erk / (erk + mek_km_erk) + erk_pp * erk_dephos;
//     dmek = mek_pp * mek_dephos - raf_kcat_mek * raf_pp * mek / (mek + raf_km_mek);
//     draf = raf_dephos * raf_pp + erk_vmax_raf * raf / (erk_km_raf * (1 + erk_pp / erk_kI_raf) + raf);

//     return true;
// }

// bool LikelihoodODE::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
// {
//     // Calculate the right-hand side of the differential equation
//     // The current timepoint and values of the dynamic variables are provided in t and y
//     // The derivative should be stored in dydt

//     Real egf_k_rafpp = parameter_values[0];
//     Real erkpp_k_raf = parameter_values[1];
//     Real rafpp_k_mekpp = parameter_values[2];
//     Real mekpp_dephos = parameter_values[3];
//     Real mekpp_k_erkpp = parameter_values[4];
//     Real dusp_k_erk = parameter_values[5];
//     Real dusp_create = parameter_values[6];
//     Real dusp_break = parameter_values[7];
//     Real egf = parameter_values[8];

//     Real erkpp = y[0];
//     Real mekpp = y[1];
//     Real rafpp = y[2];
//     Real erk = y[3];
//     Real mek = y[4];
//     Real raf = y[5];
//     Real dusp = y[6];
    
//     Real& derkpp = dydt[0];
//     Real& dmekpp = dydt[1];
//     Real& drafpp = dydt[2];
//     Real& derk = dydt[3];
//     Real& dmek = dydt[4];
//     Real& draf = dydt[5];
//     Real& ddusp = dydt[6];
    
//     derkpp = erk * mekpp * mekpp_k_erkpp - erkpp * dusp * dusp_k_erk;
//     dmekpp = mek * rafpp * rafpp_k_mekpp - mekpp * mekpp_dephos;
//     drafpp = raf * egf * egf_k_rafpp - rafpp * erkpp * erkpp_k_raf;
//     derk = erkpp * dusp * dusp_k_erk - erk * mekpp * mekpp_k_erkpp;
//     dmek = mekpp * mekpp_dephos - mek * rafpp * rafpp_k_mekpp;
//     draf = rafpp * erkpp * erkpp_k_raf - raf * egf * egf_k_rafpp;
//     ddusp = dusp_create * erkpp - dusp_break * dusp;
    
//     return true;
// }

// bool LikelihoodODE::CalculateDerivative_v6(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
// {
//     // Calculate the right-hand side of the differential equation
//     // The current timepoint and values of the dynamic variables are provided in t and y
//     // The derivative should be stored in dydt

//     Real mekpp_k_erkpp = parameter_values[0];
//     Real dusp_k_erkp = parameter_values[1];
//     Real rafpp_k_mekpp = parameter_values[2];
//     Real mekpp_dephos = parameter_values[3];
//     Real egf_k_rafpp = parameter_values[4];
//     Real erkpp_k_rafp = parameter_values[5];
//     Real mekpp_k_erkp = parameter_values[6];
//     Real dusp_k_erk = parameter_values[7];
//     Real rafpp_k_mekp = parameter_values[8];
//     Real mekp_dephos = parameter_values[9];
//     Real egf_k_rafp = parameter_values[10];
//     Real erkpp_k_raf = parameter_values[11];
//     Real egf = parameter_values[12];
//     Real erkpp_k_dusp = parameter_values[13];
//     Real dusp_break = parameter_values[14];

//     Real erkpp = y[0];
//     Real mekpp = y[1];
//     Real rafpp = y[2];
//     Real erkp = y[3];
//     Real mekp = y[4];
//     Real rafp = y[5];
//     Real erk = y[6];
//     Real mek = y[7];
//     Real raf = y[8];
//     Real dusp = y[9];
    
//     Real& derkpp = dydt[0];
//     Real& dmekpp = dydt[1];
//     Real& drafpp = dydt[2];
//     Real& derkp = dydt[3];
//     Real& dmekp = dydt[4];
//     Real& drafp = dydt[5];
//     Real& derk = dydt[6];
//     Real& dmek = dydt[7];
//     Real& draf = dydt[8];
//     Real& ddusp = dydt[9];

    
    
//     derkpp = erkp * mekpp * mekpp_k_erkpp - erkpp * dusp * dusp_k_erkp;
//     dmekpp = mekp * rafpp * rafpp_k_mekpp - mekpp * mekpp_dephos;
//     drafpp = rafp * egf * egf_k_rafpp - rafpp * erkpp * erkpp_k_rafp;
//     derkp = erkpp * dusp * dusp_k_erkp - erkp * mekpp * mekpp_k_erkpp + erk * mekpp * mekpp_k_erkp - erkp * dusp * dusp_k_erk;
//     dmekp = mekpp * mekpp_dephos - mekp * rafpp * rafpp_k_mekpp + mek * rafpp * rafpp_k_mekp - mekp * mekp_dephos;
//     drafp = rafpp * erkpp * erkpp_k_rafp - rafp * egf * egf_k_rafpp + raf * egf * egf_k_rafp - rafp * erkpp * erkpp_k_raf;
//     derk = erkp * dusp * dusp_k_erk - erk * mekpp * mekpp_k_erkp;
//     dmek = mekp * mekp_dephos - mek * rafpp * rafpp_k_mekp;
//     draf = rafp * erkpp * erkpp_k_raf - raf * egf * egf_k_rafp;
//     ddusp = erkpp * erkpp_k_dusp - dusp * dusp_break;

//     return true;
// }

bool LikelihoodODE::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
    // Calculate the right-hand side of the differential equation
    // The current timepoint and values of the dynamic variables are provided in t and y
    // The derivative should be stored in dydt
    // v7

    Real mekpp_kcat_erkpp = parameter_values[0];
    Real mekpp_km_erkpp = parameter_values[1];
    Real erkpp_dephos  = parameter_values[2];
    Real rafp_kcat_mekpp = parameter_values[3];
    Real rafp_km_mekpp = parameter_values[4];
    Real mekpp_dephos  = parameter_values[5];
    Real egf = parameter_values[7];
    Real mekpp_kcat_erkp  = parameter_values[8];
    Real mekpp_km_erkp = parameter_values[9];
    Real erkp_dephos = parameter_values[10];
    Real rafp_kcat_mekp  = parameter_values[11];
    Real rafp_km_mekp = parameter_values[12];
    Real mekp_dephos = parameter_values[13];
    Real egf_kcat_rafp  = parameter_values[14];
    Real egf_km_rafp = parameter_values[15];
    Real erkpp_kcat_raf = parameter_values[16];
    Real erkpp_km_raf  = parameter_values[17];

    Real erkpp = y[0];
    Real mekpp = y[1];
    Real erkp = y[2];
    Real mekp = y[3];
    Real rafp = y[4];
    Real erk = y[5];
    Real mek = y[6];
    Real raf = y[7];
    
    Real& derkpp = dydt[0];
    Real& dmekpp = dydt[1];
    Real& derkp = dydt[2];
    Real& dmekp = dydt[3];
    Real& drafp = dydt[4];
    Real& derk = dydt[5];
    Real& dmek = dydt[6];
    Real& draf = dydt[7];

    derkpp = mekpp_kcat_erkpp * mekpp * erkp / (erkp + mekpp_km_erkpp) - erkpp * erkpp_dephos;
    dmekpp = rafp_kcat_mekpp * rafp * mekp / (rafp_km_mekpp + mekp) - mekpp * mekpp_dephos;
    derkp = erkpp * erkpp_dephos - mekpp_kcat_erkpp * mekpp * erkp / (erkp + mekpp_km_erkpp) + mekpp_kcat_erkp * mekpp * erk / (mekpp_km_erkp + erk) - erkp * erkp_dephos;
    dmekp = mekpp * mekpp_dephos - rafp_kcat_mekpp * rafp * mekp / (rafp_km_mekpp + mekp) + rafp_kcat_mekp * rafp * mek / (rafp_km_mekp + mek) - mekp * mekp_dephos;
    drafp = egf_kcat_rafp * egf * raf / (egf_km_rafp + raf) - erkpp_kcat_raf * erkpp * rafp / (erkpp_km_raf + rafp);
    derk = erkp * erkp_dephos - mekpp_kcat_erkp * mekpp * erk / (mekpp_km_erkp + erk);
    dmek = mekp * mekp_dephos - rafp_kcat_mekp * rafp * mek / (rafp_km_mekp + mek);
    draf = erkpp_kcat_raf * erkpp * rafp / (erkpp_km_raf + rafp) - egf_kcat_rafp * egf * raf / (egf_km_rafp + raf);

    return true;
}


