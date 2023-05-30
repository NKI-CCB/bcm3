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
    size_t num_dynamic_variables = 6;
    size_t num_inference_variables = 16;

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
    parser.Parse(cwd.string(), ",", false);

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

    // // The initial conditions can be dependent on the parameters
    // OdeVectorReal initial_conditions(4);
    // initial_conditions(0) = parameter_values(9);
    // initial_conditions(1) = parameter_values(10);
    // initial_conditions(2) = parameter_values(11);
    // initial_conditions(3) = parameter_values(12);

    // The initial conditions can be dependent on the parameters
    OdeVectorReal initial_conditions(6);
    initial_conditions(0) = parameter_values(10);
    initial_conditions(1) = parameter_values(11);
    initial_conditions(2) = parameter_values(12);
    initial_conditions(3) = parameter_values(13);
    initial_conditions(4) = parameter_values(14);
    initial_conditions(5) = parameter_values(15);

    boost::filesystem::path cwd = boost::filesystem::current_path() / "normalized_oscillations.csv";
    
    bcm3::CSVParser parser;
    parser.Parse("/Users/huubvdent/Documents/Internship/BCM_RUNS/MAPK_ERK_v4/data_10_fit/normalized_oscillations.csv", ",", false);

    // Integrate the ODE system and calculate the likelihood based on the solution
    // In this example we compare the first dynamic variable to a cosine
    if (solver->Simulate(initial_conditions.data(), timepoints, simulated_trajectories)) {
        logp = 0.0;
        for (size_t i = 0; i < timepoints.size(); i++) {
            // Real cosvalue = 100.0 * cos(timepoints(i) / 1140) + 150.0;
            Real datavalue = parser.GetEntry(10, i);
            logp += bcm3::LogPdfTnu3(datavalue, simulated_trajectories(0, i), 20.0);
        }
    } else {
        logp = -std::numeric_limits<Real>::infinity();
    }

	return true;
}

// bool LikelihoodODE::CalculateDerivative_v3(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
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

bool LikelihoodODE::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
    // Calculate the right-hand side of the differential equation
    // The current timepoint and values of the dynamic variables are provided in t and y
    // The derivative should be stored in dydt

    Real erk_dephos = parameter_values[0];
    Real erk_kI_raf = parameter_values[1];
    Real erk_km_raf = parameter_values[2];
    Real erk_vmax_raf = parameter_values[3];
    Real mek_dephos = parameter_values[4];
    Real mek_kcat_erk = parameter_values[5];
    Real mek_km_erk = parameter_values[6];
    Real raf_dephos = parameter_values[7];
    Real raf_kcat_mek = parameter_values[8];
    Real raf_km_mek = parameter_values[9];

    Real erk_pp = y[0];
    Real mek_pp = y[1];
    Real raf_pp = y[2];
    Real erk = y[3];
    Real mek = y[4];
    Real raf = y[5];

    Real& derk_pp = dydt[0];
    Real& dmek_pp = dydt[1];
    Real& draf_pp = dydt[2];
    Real& derk = dydt[3];
    Real& dmek = dydt[4];
    Real& draf = dydt[5];
    
    derk_pp = mek_kcat_erk * mek_pp * erk / (erk + mek_km_erk) - erk_pp * erk_dephos;
    dmek_pp = - mek_pp * mek_dephos + raf_kcat_mek * raf_pp * mek /(mek + raf_km_mek);
    draf_pp = raf_dephos * raf_pp - erk_vmax_raf * raf / (erk_km_raf * (1 + erk_pp / erk_kI_raf) + raf);
    derk = mek_kcat_erk * mek_pp * erk / (erk + mek_km_erk) + erk_pp * erk_dephos;
    dmek = mek_pp * mek_dephos - raf_kcat_mek * raf_pp * mek / (mek + raf_km_mek);
    draf = raf_dephos * raf_pp + erk_vmax_raf * raf / (erk_km_raf * (1 + erk_pp / erk_kI_raf) + raf);

    return true;
}
