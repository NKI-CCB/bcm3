#include "Utils.h"
#include "TestLikelihoodTruncatedT.h"
#include "ProbabilityDistributions.h"
#include "VectorUtils.h"
#include "mvt.h"

TestLikelihoodTruncatedT::TestLikelihoodTruncatedT(size_t sampling_threads, size_t evaluation_threads)
	: dimensions(0)
	, num_clusters(0)
{
}

TestLikelihoodTruncatedT::~TestLikelihoodTruncatedT()
{
}

bool TestLikelihoodTruncatedT::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	try {
		dimensions = likelihood_node.get<size_t>("<xmlattr>.dimensions");
		num_clusters = likelihood_node.get<size_t>("<xmlattr>.num_clusters");

		if (varset->GetNumVariables() != dimensions) {
			LOGERROR("Incorrect number of variables in prior for samples a %z-dimensional space", dimensions);
			return false;
		}

		mus.resize(num_clusters, VectorReal::Constant(dimensions, std::numeric_limits<Real>::quiet_NaN()));
		sigmas.resize(num_clusters, MatrixReal::Constant(dimensions, dimensions, std::numeric_limits<Real>::quiet_NaN()));
		for (size_t i = 0; i < num_clusters; i++) {
			std::string mu = likelihood_node.get<std::string>(std::string("<xmlattr>.mu") + std::to_string(i+1));
			if (!bcm3::ParseVectorFromString(mu, mus[i])) {
				return false;
			}
			if (mus[i].size() != dimensions) {
				LOGERROR("Inconsistent dimension for mu%z", i);
				return false;
			}

			std::string sigma = likelihood_node.get<std::string>(std::string("<xmlattr>.sigma") + std::to_string(i + 1));
			if (!bcm3::ParseMatrixFromString(sigma, sigmas[i])) {
				return false;
			}
			if (sigmas[i].rows() != dimensions || sigmas[i].cols() != dimensions) {
				LOGERROR("Inconsistent dimension for sigma%z", i);
				return false;
			}
		}

		std::string nustr = likelihood_node.get<std::string>(std::string("<xmlattr>.nus"));
		if (!bcm3::ParseVectorFromString(nustr, nus)) {
			return false;
		}
		if (nus.size() != num_clusters) {
			LOGERROR("Inconsistent number of nus");
			return false;
		}

		std::string weightstr = likelihood_node.get<std::string>(std::string("<xmlattr>.weights"));
		if (!bcm3::ParseVectorFromString(weightstr, weights)) {
			return false;
		}
		if (weights.size() != num_clusters) {
			LOGERROR("Inconsistent number of weights");
			return false;
		}
		weights.array() /= weights.sum();
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return true;
}

bool TestLikelihoodTruncatedT::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = -std::numeric_limits<Real>::infinity();

	for (size_t i = 0; i < num_clusters; i++) {
		Real thislogp = bcm3::dmvt(values, mus[i], sigmas[i], nus[i], true);
		logp = bcm3::logsum(logp, thislogp + log(weights[i]));
	}

	return true;
}
