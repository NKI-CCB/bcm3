#include "Utils.h"
#include "PharmacoLikelihoodSingle.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include <boost/property_tree/xml_parser.hpp>

PharmacoLikelihoodSingle::PharmacoLikelihoodSingle(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, evaluation_threads(evaluation_threads)
	, additive_sd_ix(std::numeric_limits<size_t>::max())
	, proportional_sd_ix(std::numeric_limits<size_t>::max())
{
}

PharmacoLikelihoodSingle::~PharmacoLikelihoodSingle()
{
}

bool PharmacoLikelihoodSingle::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	bool result = true;;

#if 0
	try {
		// Load all experiments
		BOOST_FOREACH(const boost::property_tree::ptree::value_type& data, likelihood_node.get_child("")) {
			if (data.first == "experiment") {
				std::unique_ptr<Experiment> experiment = Experiment::Create(data.second, varset, vm, rng, evaluation_threads, store_simulation);
				if (!experiment) {
					return false;
				}
				experiments.push_back(std::move(experiment));
			}
		}
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}
#endif

	return result;
}

bool PharmacoLikelihoodSingle::PostInitialize()
{
	additive_sd_ix = varset->GetVariableIndex("additive_error_standard_deviation");
	proportional_sd_ix = varset->GetVariableIndex("proportional_error_standard_deviation");
	if (additive_sd_ix == std::numeric_limits<size_t>::max() && proportional_sd_ix == std::numeric_limits<size_t>::max()) {
		LOGERROR("Neither \"additive_error_standard_deviation\" nor \"proportional_error_standard_deviation\" has been specified in the prior; at least one of these variables should be included.");
		return false;
	}

	return true;
}

bool PharmacoLikelihoodSingle::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	Real additive_sd = 0.0;
	if (additive_sd_ix != std::numeric_limits<size_t>::max()) {
		additive_sd = varset->TransformVariable(additive_sd_ix, values(additive_sd_ix));
	}
	Real proportional_sd = 0.0;
	if (proportional_sd_ix != std::numeric_limits<size_t>::max()) {
		varset->TransformVariable(proportional_sd_ix, values(proportional_sd_ix));
	}

	Real conversion = 1.0;

	if (model.Solve(treatment_timepoints, treatment_doses, observation_timepoints, simulated_data)) {
		for (ptrdiff_t i = 0; i < observation_timepoints.size(); i++) {
			Real x = conversion * simulated_data(i);
			Real y = observed_data(i);

			if (!std::isnan(y)) {
				logp += bcm3::LogPdfTnu4(x, y, additive_sd + proportional_sd * std::max(x, 0.0));
			}
		}
	} else {
		logp = -std::numeric_limits<Real>::infinity();
	}

	return true;
}

void PharmacoLikelihoodSingle::AddOptionsDescription(boost::program_options::options_description& pod)
{
}
