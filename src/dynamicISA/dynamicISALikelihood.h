#pragma once

#include "Likelihood.h"

class dynamicISAExperiment;
class SignalingModel;

class dynamicISALikelihood : public bcm3::Likelihood
{
public:
	dynamicISALikelihood(size_t numthreads, size_t evaluation_threads);
	~dynamicISALikelihood();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node);
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	size_t GetNumSignalingMolecules() const;
	const dynamicISAExperiment* GetExperiment(const std::string& experiment) const;

private:
	std::shared_ptr<SignalingModel> model;
	std::shared_ptr<const bcm3::VariableSet> varset;
	std::vector< std::unique_ptr<dynamicISAExperiment> > experiments;

	VectorReal transformed_variables;
};
