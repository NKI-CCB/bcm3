#pragma once

#include "ODESolver.h"
#include "Likelihood.h"

class LikelihoodODE : public bcm3::Likelihood
{
public:
	LikelihoodODE(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodODE();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	const OdeMatrixReal& GetSimulatedTrajectories() const { return simulated_trajectories; }

private:
	bool CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);

	std::shared_ptr<const bcm3::VariableSet> varset;
	std::shared_ptr<ODESolver> solver;
	OdeVectorReal timepoints;
	VectorReal parameter_values;
	OdeMatrixReal simulated_trajectories;
};
