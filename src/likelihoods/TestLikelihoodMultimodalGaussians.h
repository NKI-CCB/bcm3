#pragma once

#include "Likelihood.h"

class TestLikelihoodMultimodalGaussians : public bcm3::Likelihood
{
public:
	TestLikelihoodMultimodalGaussians(size_t sampling_threads, size_t evaluation_threads);
	~TestLikelihoodMultimodalGaussians();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	std::shared_ptr<const bcm3::VariableSet> varset;
	size_t dimension;

	std::vector<VectorReal> means;
	std::vector<MatrixReal> covariances;
};
