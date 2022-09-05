#pragma once

#include "VariableSet.h"

#include <boost/program_options.hpp>

namespace bcm3 {

class Likelihood
{
public:
	virtual ~Likelihood();

	bool SetLearningRate(Real learning_rate);
	inline Real GetLearningRate() const { return learning_rate; }

	virtual bool Initialize(std::shared_ptr<const VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool AddNonSampledParameters(const std::vector<std::string>& variable_names);
	virtual void SetNonSampledParameters(const VectorReal& values);
	virtual bool PostInitialize();
	virtual bool IsReentrant() = 0;
	virtual void OutputEvaluationStatistics(const std::string& path) const {}

	//! Evaluate the log likelihood for the given values of the variables.
	//! @param threadix An index indicating which thread is calling the function.
	//! @param values The values of the variables, corresponding to the VariableSet provided to the sampler.
	//! @param logp The log likelihood value should be written into this parameter.
	//! @return True if successful, false if a non-recoverable error occured (this will cause the sampler to stop).
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp) = 0;

protected:
	Likelihood();

	Real learning_rate;
};

}
