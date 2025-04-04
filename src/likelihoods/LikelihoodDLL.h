#pragma once

#include "Likelihood.h"

class LikelihoodDLL : public bcm3::Likelihood
{
public:
	LikelihoodDLL(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodDLL();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize();
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	typedef bool (*initialize_fn)(size_t num_variables, const char* const* variable_names);
	typedef bool (*likelihood_fn)(size_t num_variables, const double* values, const char* const* variable_names, double* log_p);

#if PLATFORM_WINDOWS
	HMODULE likelihood_dll;
#else
	void* likelihood_dll;
#endif

	initialize_fn initialize;
	likelihood_fn likelihood;
	size_t num_variables;
	char** variable_names;
};
