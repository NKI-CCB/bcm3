#pragma once

#include "Likelihood.h"

class LikelihoodMitosisTimeEstimation : public bcm3::Likelihood
{
public:
	LikelihoodMitosisTimeEstimation(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodMitosisTimeEstimation();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	std::shared_ptr<const bcm3::VariableSet> varset;

	VectorReal timepoints;
	MatrixReal observed_trajectories;
	MatrixReal sobol_sequence_values;
	Real misspecification_time;

	// Runtime variables
	MatrixReal simulated_trajectories;
	MatrixReal likelihoods;
};
