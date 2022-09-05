#include "Utils.h"
#include "Likelihood.h"

namespace bcm3 {

Likelihood::Likelihood()
	: learning_rate(1.0)
{
}

Likelihood::~Likelihood()
{
}

bool Likelihood::SetLearningRate(Real learning_rate)
{
	if (learning_rate < 0.0 || learning_rate > 1.0) {
		LOGERROR("Learning rate must be >= 0 and <= 1.0");
		return false;
	}
	this->learning_rate = learning_rate;
	return true;
}

bool Likelihood::Initialize(std::shared_ptr<const VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	return true;
}

bool Likelihood::AddNonSampledParameters(const std::vector<std::string>& variable_names)
{
	return true;
}

void Likelihood::SetNonSampledParameters(const VectorReal& values)
{
}

bool Likelihood::PostInitialize()
{
	return true;
}

}
