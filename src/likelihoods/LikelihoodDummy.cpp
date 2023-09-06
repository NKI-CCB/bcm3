#include "Utils.h"
#include "LikelihoodDummy.h"

LikelihoodDummy::LikelihoodDummy(size_t sampling_threads, size_t evaluation_threads)
{
}

LikelihoodDummy::~LikelihoodDummy()
{
}

bool LikelihoodDummy::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	return true;
}

bool LikelihoodDummy::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;
	return true;
}
