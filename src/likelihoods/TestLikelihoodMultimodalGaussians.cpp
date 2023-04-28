#include "Utils.h"
#include "TestLikelihoodMultimodalGaussians.h"
#include "ProbabilityDistributions.h"
#include "mvn.h"

TestLikelihoodMultimodalGaussians::TestLikelihoodMultimodalGaussians(size_t sampling_threads, size_t evaluation_threads)
{
}

TestLikelihoodMultimodalGaussians::~TestLikelihoodMultimodalGaussians()
{
}

bool TestLikelihoodMultimodalGaussians::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;
	
	if (varset->GetNumVariables() != 2) {
		LOGERROR("Inconsistent prior and likelihood (%zu variables and %zu dimensions)", varset->GetNumVariables(), 2);
		return false;
	}

	means.resize(2, VectorReal::Zero(2));
	means[0] << -5, -5;
	means[1] << 5, 5;

	covariances.resize(2, MatrixReal::Identity(2, 2));
	covariances[0] << 1, -0.5, -0.5, 0.5;
	covariances[1] << 2, 1, 1, 1; 

	return true;
}

bool TestLikelihoodMultimodalGaussians::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	Real logp1 = log(0.5) + bcm3::dmvnormal(values, means[0], covariances[0], true);
	Real logp2 = log(0.5) + bcm3::dmvnormal(values, means[1], covariances[1], true);

	logp = bcm3::logsum(logp1, logp2);

	return true;
}
