#pragma once

#include "Likelihood.h"
#include "VariableSet.h"

#include <boost/program_options.hpp>

namespace bcm3 {

class LikelihoodFactory
{
public:
	static std::shared_ptr<bcm3::Likelihood> CreateLikelihood(std::string likelihood_xml_fn, std::shared_ptr<const bcm3::VariableSet> varset, const boost::program_options::variables_map& vm, size_t sampling_threads, size_t evaluation_threads);

	static void AddOptionsDescription(boost::program_options::options_description& pod);

private:
	LikelihoodFactory() {}
	~LikelihoodFactory() {}
};

}
