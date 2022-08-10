#include "Utils.h"
#include "TestLikelihoodCircular.h"
#include "ProbabilityDistributions.h"

TestLikelihoodCircular::TestLikelihoodCircular(size_t sampling_threads, size_t evaluation_threads)
{
}

TestLikelihoodCircular::~TestLikelihoodCircular()
{
}

bool TestLikelihoodCircular::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node)
{
	this->varset = varset;

	try {
		boost::property_tree::ptree model_node = likelihood_node.get_child("model");
		std::string model_fn = model_node.get<std::string>("<xmlattr>.file");
		//dim = model_node.get<size_t>("<xmlattr>.dimension");
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return true;
}

bool TestLikelihoodCircular::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	return true;
}
