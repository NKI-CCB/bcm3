#include "Utils.h"
#include "TestLikelihoodBanana.h"
#include "ProbabilityDistributions.h"

TestLikelihoodBanana::TestLikelihoodBanana(size_t sampling_threads, size_t evaluation_threads)
{
}

TestLikelihoodBanana::~TestLikelihoodBanana()
{
}

bool TestLikelihoodBanana::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	try {
		dim = likelihood_node.get<size_t>("<xmlattr>.dimension");
		if (dim != varset->GetNumVariables()) {
			LOGERROR("Dimension %zd does not match number of variables in prior (%zd)", dim, varset->GetNumVariables());
			return false;
		}
		if (dim < 2) {
			LOGERROR("Dimension is %zd but should be at least 2", dim);
			return false;
		}

		sd1 = likelihood_node.get<Real>("<xmlattr>.sd1");
		sd2 = likelihood_node.get<Real>("<xmlattr>.sd2");
		if (sd1 <= 0.0 || sd2 <= 0.0) {
			LOGERROR("Standard deviations should be greater than 0");
			return false;
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return true;
}

bool TestLikelihoodBanana::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	Real p = 1.0;
	for (size_t i = 0; i < dim-1; i++) {
		p = p * bcm3::PdfNormal(values(i), 0, sd1);
	}
	Real y = values(0);
	for (size_t i = 1; i < dim - 1; i++) {
		y += values(i);
	}
	p *= bcm3::PdfNormal(values(dim-1), y+3*y+bcm3::square(1-y), sd2);
	logp = log(p);
	return true;
}
