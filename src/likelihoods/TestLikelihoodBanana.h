#pragma once

#include "Likelihood.h"

class TestLikelihoodBanana : public bcm3::Likelihood
{
public:
	TestLikelihoodBanana(size_t sampling_threads, size_t evaluation_threads);
	~TestLikelihoodBanana();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	std::shared_ptr<const bcm3::VariableSet> varset;
	size_t dim;
	Real sd1;
	Real sd2;
};
