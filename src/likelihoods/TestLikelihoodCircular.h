#pragma once

#include "Likelihood.h"

class TestLikelihoodCircular : public bcm3::Likelihood
{
public:
	TestLikelihoodCircular(size_t sampling_threads, size_t evaluation_threads);
	~TestLikelihoodCircular();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node);
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	std::shared_ptr<const bcm3::VariableSet> varset;
};
