#include "Utils.h"
#include "LikelihoodODE.h"
#include "ODESolverCVODE.h"
#include "ProbabilityDistributions.h"

LikelihoodODE::LikelihoodODE(size_t sampling_threads, size_t evaluation_threads)
{
}

LikelihoodODE::~LikelihoodODE()
{
}

bool LikelihoodODE::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
    size_t num_dynamic_variables = 4;
    size_t num_inference_variables = 13;

    this->varset = varset;
    if (varset->GetNumVariables() != num_inference_variables) {
        LOGERROR("Incorrect number of parameters");
        return false;
    }

    solver = std::make_shared<ODESolverCVODE>();

	ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodODE::CalculateDerivative, this, _1, _2, _3, _4);
	solver->SetDerivativeFunction(derivative);
	solver->Initialize(num_dynamic_variables, NULL);
	solver->SetTolerance(1e-8, 1e-8);

    // Allocate parameter vector to number of parameters that are used
    parameter_values.resize(num_inference_variables);

    // Timepoints at which the ODE solver should return values of the dynamic variables
    Real max_time = 1000.0;
    size_t num_timepoints = 100;
    timepoints.resize(num_timepoints);
    for (size_t i = 0; i < num_timepoints; i++) {
        timepoints(i) = 0 + max_time * (i / (Real)(num_timepoints - 1));
    }

	return true;
}

bool LikelihoodODE::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
    // Apply variable transformations if those were specified in the prior
    for (ptrdiff_t i = 0; i < values.size(); i++) {
        parameter_values(i) = varset->TransformVariable(i, values(i));
    }

    // The initial conditions can be dependent on the parameters
    OdeVectorReal initial_conditions(solver->GetNumVariables());
    initial_conditions(0) = parameter_values(9);
    initial_conditions(1) = parameter_values(10);
    initial_conditions(2) = parameter_values(11);
    initial_conditions(3) = parameter_values(12);

    // Integrate the ODE system and calculate the likelihood based on the solution
    // In this example we compare the first dynamic variable to a cosine
    if (solver->Simulate(initial_conditions.data(), timepoints, simulated_trajectories)) {
        logp = 0.0;
        for (size_t i = 0; i < timepoints.size(); i++) {
            Real cosvalue = 100.0 * cos(timepoints(i) / 2300.0) + 300.0;
            logp += bcm3::LogPdfTnu3(cosvalue, simulated_trajectories(0, i), 10.0);
        }
    } else {
        logp = -std::numeric_limits<Real>::infinity();
    }

	return true;
}

bool LikelihoodODE::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
    // Calculate the right-hand side of the differential equation
    // The current timepoint and values of the dynamic variables are provided in t and y
    // The derivative should be stored in dydt

    return true;
}
