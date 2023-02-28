#pragma once

#include "Likelihood.h"

class TestLikelihoodTruncatedT : public bcm3::Likelihood
{
public:
	TestLikelihoodTruncatedT(size_t sampling_threads, size_t evaluation_threads);
	~TestLikelihoodTruncatedT();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	std::shared_ptr<const bcm3::VariableSet> varset;

	size_t dimensions;
	size_t num_clusters;
	VectorReal nus;
	VectorReal weights;
	std::vector<VectorReal> mus;
	std::vector<MatrixReal> sigmas;
};
