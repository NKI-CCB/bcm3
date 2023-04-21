#include "Utils.h"
#include "TestLikelihoodCircular.h"
#include "ProbabilityDistributions.h"

TestLikelihoodCircular::TestLikelihoodCircular(size_t sampling_threads, size_t evaluation_threads)
{
}

TestLikelihoodCircular::~TestLikelihoodCircular()
{
}

bool TestLikelihoodCircular::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	Real offset = std::numeric_limits<Real>::quiet_NaN();

	try {
		dimension = likelihood_node.get<size_t>("<xmlattr>.dimension");
		r = likelihood_node.get<Real>("<xmlattr>.radius", 2.0);
		offset = likelihood_node.get<Real>("<xmlattr>.offset", 3.5);
		w = likelihood_node.get<Real>("<xmlattr>.width", 0.1);
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	if (varset->GetNumVariables() != dimension) {
		LOGERROR("Inconsistent prior and likelihood (%zu variables and %zu dimensions)", varset->GetNumVariables(), dimension);
		return false;
	}

	mu1.setZero(dimension);
	mu2.setZero(dimension);
	mu1(0) = -offset;
	mu2(0) = offset;

	return true;
}

bool TestLikelihoodCircular::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	Real x1 = (values - mu1).norm();
	Real x2 = (values - mu2).norm();

	Real logp1 = bcm3::LogPdfNormal(x1, r, w);
	Real logp2 = bcm3::LogPdfNormal(x2, r, w);

	logp = bcm3::logsum(logp1, logp2);

	return true;
}
