#pragma once

#include "Likelihood.h"
#include "VariableSet.h"

class LikelihoodFactory
{
public:
	static std::shared_ptr<bcm3::Likelihood> CreateLikelihood(std::string likelihood_xml_fn, std::shared_ptr<const bcm3::VariableSet> varset, size_t sampling_threads, size_t evaluation_threads);

private:
	LikelihoodFactory() {}
	~LikelihoodFactory() {}
};
