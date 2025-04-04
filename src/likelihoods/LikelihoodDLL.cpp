#include "Utils.h"
#include "LikelihoodDLL.h"

LikelihoodDLL::LikelihoodDLL(size_t sampling_threads, size_t evaluation_threads)
	: likelihood_dll(nullptr)
	, likelihood(nullptr)
	, num_variables(0)
	, variable_names(nullptr)
{
}

LikelihoodDLL::~LikelihoodDLL()
{
	if (variable_names != nullptr) {
		for (ptrdiff_t i = 0; i < num_variables; i++) {
			delete[] variable_names[i];
		}
		delete[] variable_names;
	}

	if (likelihood_dll) {
#if PLATFORM_WINDOWS
		FreeLibrary(likelihood_dll);
#else
		dlclose(likelihood_dll);
#endif
	}
}

bool LikelihoodDLL::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	num_variables = varset->GetNumVariables();
	variable_names = new char* [num_variables];
	for (size_t i = 0; i < num_variables; i++) {
		variable_names[i] = new char[varset->GetVariableName(i).size() + 1];
		strcpy(variable_names[i], varset->GetVariableName(i).c_str());
	}

	std::string dll_filename;
	bool include_build_dir;
	try {
		dll_filename = likelihood_node.get<std::string>("<xmlattr>.dll_filename_base");
		include_build_dir = likelihood_node.get<bool>("<xmlattr>.include_build_dir", true);
		
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

#if PLATFORM_WINDOWS
	dll_filename += ".dll";
	if (include_build_dir) {
#if BUILD_DEBUG
		dll_filename = "build/Debug/" + dll_filename;
#else
		dll_filename = "build/Release/" + dll_filename;
#endif
	}
	likelihood_dll = LoadLibrary(dll_filename.c_str());
	if (likelihood_dll) {
		initialize = (initialize_fn)GetProcAddress(likelihood_dll, "initialize_likelihood");
		likelihood = (likelihood_fn)GetProcAddress(likelihood_dll, "evaluate_log_probability");
#else
	dll_filename = "build/" + dll_filename + ".so";
	if (add_build_type) {
	likelihood_dll = dlopen(dll_filename.c_str(), RTLD_NOW);
	if (likelihood_dll) {
		initialize = (initialize_fn)dlsym(likelihood_dll, "initialize_likelihood");
		likelihood = (likelihood_fn)dlsym(likelihood_dll, "evaluate_log_probability");
#endif

		if (!likelihood) {
			LOGERROR("Unable to find evaluate_log_probability function in dll \"%s\".", dll_filename.c_str());
			return false;
		}
	} else {
		LOGERROR("Can't find dll for likelihood function \"%s\".", dll_filename.c_str());
		return false;
	}

	return true;
}

bool LikelihoodDLL::PostInitialize()
{
	if (initialize) {
		if (!initialize(num_variables, variable_names)) {
			LOGERROR("DLL initialize function returning false; halting inference.");
			return false;
		}
	}
	return true;
}

bool LikelihoodDLL::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	double logp_double = std::numeric_limits<double>::quiet_NaN();
	bool result = likelihood(num_variables, values.data(), variable_names, &logp_double);
	if (result) {
		logp = logp_double;
	} else {
		return false;
	}

	if (std::isnan(logp)) {
		return false;
	} else {
		return true;
	}
}
