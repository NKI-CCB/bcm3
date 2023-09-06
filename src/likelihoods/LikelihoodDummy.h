#pragma once

#include "Likelihood.h"

class LikelihoodDummy : public bcm3::Likelihood
{
public:
	LikelihoodDummy(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodDummy();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
};
