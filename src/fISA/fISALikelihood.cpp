#include "Utils.h"
#include "fISALikelihood.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"
#include "SignalingNetwork.h"
#include "fISAExperiment.h"

#include <boost/property_tree/xml_parser.hpp>

fISALikelihood::fISALikelihood(size_t numthreads, size_t evaluation_threads)
	: numthreads(numthreads)
	, evaluation_threads(evaluation_threads)
	, cell_line_evaluation_task_manager(evaluation_threads)
{
}

fISALikelihood::~fISALikelihood()
{
}

bool fISALikelihood::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	bool result = true;

	this->VarSet = varset;
	
	// Load the model
	boost::property_tree::ptree model = likelihood_node.get_child("model");
	std::string modelfn = model.get<std::string>("<xmlattr>.file");
	std::string activation_limit_str = model.get<std::string>("<xmlattr>.activation_limit", "logistic");
	size_t multiroot_solves = model.get<size_t>("<xmlattr>.multiroot_solves", 10);

	SignalingNetwork::ActivationLimitType activation_limit;
	if (activation_limit_str == "minmax") {
		activation_limit = SignalingNetwork::ActivationLimit_MinMax;
	} else if (activation_limit_str == "logistic") {
		activation_limit = SignalingNetwork::ActivationLimit_Logistic;
	} else if (activation_limit_str == "logistic_or") {
		activation_limit = SignalingNetwork::ActivationLimit_LogisticOr;
	} else {
		LOGERROR("Unknown activation limit \"%s\"", activation_limit_str.c_str());
		return false;
	}

	network = std::make_unique<SignalingNetwork>(evaluation_threads);
	if (!network->Initialize(modelfn, VarSet, activation_limit, multiroot_solves)) {
		return false;
	}
		
	// Load all experiments
	BOOST_FOREACH(const boost::property_tree::ptree::value_type& data, likelihood_node.get_child("")) {
		if (data.first == "experiment") {
			std::unique_ptr<fISAExperiment> experiment = fISAExperiment::Create(data.second, network.get(), VarSet, experiments, evaluation_threads, &cell_line_evaluation_task_manager);
			if (!experiment) {
				return false;
			}
			experiments.push_back(std::move(experiment));
		}
	}

	return result;
}

#if 0
void fISALikelihood::SetBootstrap(unsigned long rngseed)
{
	for (size_t i = 0; i < experiments.size(); i++) {
		experiments[i]->SetBootstrap(rngseed);
	}
}

void fISALikelihood::SetSafeBayesIx(size_t ix)
{
	for (size_t i = 0; i < experiments.size(); i++) {
		experiments[i]->SetSafeBayesIx(ix);
	}
}

void fISALikelihood::AddSafeBayesCvIx(size_t ix)
{
	for (size_t i = 0; i < experiments.size(); i++) {
		experiments[i]->AddSafeBayesCvIx(ix);
	}
}
#endif

bool fISALikelihood::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	for (size_t i = 0; i < experiments.size(); i++) {
		experiments[i]->StartEvaluateLogProbability(values, experiments);
	}

	for (size_t i = 0; i < experiments.size(); i++) {
		Real experiment_logp = 0.0;
		bool result = experiments[i]->FinishEvaluateLogProbability(values, experiment_logp, experiments);
		if (result) {
			logp += experiment_logp;
		} else {
			return false;
		}
	}

	return true;
}

size_t fISALikelihood::GetNumSignalingMolecules() const
{
	return network->GetMoleculeCount();
}

const fISAExperiment* fISALikelihood::GetExperiment(const std::string& experiment) const
{
	for (size_t i = 0; i < experiments.size(); i++) {
		if (experiments[i]->GetName() == experiment) {
			return experiments[i].get();
		}
	}
	LOGERROR("Experiment \"%s\" doesn't exist", experiment.c_str());
	return nullptr;
}
