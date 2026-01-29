#include "Utils.h"
#include "LikelihoodDummy.h"
#include "ProbabilityDistributions.h"

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
#if 0
	Real x[10] = {
		2.4187462, 2.2157475, 1.0286695, 1.8008323, 2.6945451, 2.7897943, 1.1919238, 2.0363240, 0.7733273, 2.1435216
	};

	logp = 0.0;
	for (int i = 0; i < 10; i++) {
		logp += bcm3::LogPdfNormal(x[i], values[0], values[1]);
	}
#else
	logp = bcm3::LogPdfTnu4(values[0], 0, 1);
#endif
	return true;
}
