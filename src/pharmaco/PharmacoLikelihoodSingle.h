#pragma once

#include "Likelihood.h"
#include "PharmacokineticModel.h"

#include <boost/program_options.hpp>

class PharmacoLikelihoodSingle : public bcm3::Likelihood
{
public:
	PharmacoLikelihoodSingle(size_t sampling_threads, size_t evaluation_threads);
	~PharmacoLikelihoodSingle();
	
	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize();
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	static void AddOptionsDescription(boost::program_options::options_description& pod);

private:
	// Static variables
	size_t sampling_threads;
	size_t evaluation_threads;
	std::shared_ptr<const bcm3::VariableSet> varset;

	size_t additive_sd_ix;
	size_t proportional_sd_ix;

	PharmacokineticModel model;

	VectorReal treatment_timepoints;
	VectorReal treatment_doses;
	VectorReal observation_timepoints;
	VectorReal observed_data;
	VectorReal simulated_data;
};
