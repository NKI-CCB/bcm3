#include "Utils.h"
#include "CellPopulationLikelihood.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include <boost/property_tree/xml_parser.hpp>

CellPopulationLikelihood::CellPopulationLikelihood(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, evaluation_threads(evaluation_threads)
{
}

CellPopulationLikelihood::~CellPopulationLikelihood()
{
}

bool CellPopulationLikelihood::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;
	transformed_variables.setConstant(varset->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());

	bool result = true;
	try {
		// Load all experiments
		BOOST_FOREACH(const boost::property_tree::ptree::value_type& data, likelihood_node.get_child("")) {
			if (data.first == "experiment") {
				std::unique_ptr<Experiment> experiment = Experiment::Create(data.second, varset, vm, rng, evaluation_threads);
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

	experiment_likelihood = VectorReal::Constant(experiments.size(), std::numeric_limits<Real>::quiet_NaN());

	return result;
}

bool CellPopulationLikelihood::AddNonSampledParameters(const std::vector<std::string>& variable_names)
{
	for (size_t i = 0; i < experiments.size(); i++) {
		if (!experiments[i]->AddNonSampledParameters(variable_names)) {
			return false;
		}
	}
	return true;
}

void CellPopulationLikelihood::SetNonSampledParameters(const VectorReal& values)
{
	for (size_t i = 0; i < experiments.size(); i++) {
		experiments[i]->SetNonSampledParameters(values);
	}
}

bool CellPopulationLikelihood::PostInitialize()
{
	for (size_t i = 0; i < experiments.size(); i++) {
		if (!experiments[i]->PostInitialize()) {
			return false;
		}
	}
	return true;
}

void CellPopulationLikelihood::OutputEvaluationStatistics(const std::string& path) const
{
#if 0
	for (size_t i = 0; i < experiments.size(); i++) {
		experiments[i]->DumpCVodeStatistics(path);
	}
#endif
}

bool CellPopulationLikelihood::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	for (size_t i = 0; i < varset->GetNumVariables(); i++) {
		transformed_variables(i) = varset->TransformVariable(i, values[i]);
	}

	for (size_t i = 0; i < experiments.size(); i++) {
		Real experiment_logp;
		if (!experiments[i]->EvaluateLogProbability(threadix, values, transformed_variables, experiment_logp)) {
			logp = -std::numeric_limits<Real>::infinity();
			return false;
		}
		experiment_likelihood[i] = experiment_logp;
		logp += experiment_logp;
	}

	return true;
}

const Experiment* CellPopulationLikelihood::GetExperiment(const std::string& experiment)
{
	for (size_t i = 0; i < experiments.size(); i++) {
		if (experiments[i]->GetName() == experiment) {
			return experiments[i].get();
		}
	}
	return NULL;
}

void CellPopulationLikelihood::AddOptionsDescription(boost::program_options::options_description& pod)
{
	pod.add_options()
		("cellpop.use_only_cell_ix", boost::program_options::value<std::string>()->default_value("-1"), "Use only a specific cell in the data likelihood; -1 to disable.")
	;
}
