#include "Utils.h"
#include "dynamicISALikelihood.h"
#include "dynamicISAExperiment.h"
#include "SignalingModel.h"

dynamicISALikelihood::dynamicISALikelihood(size_t sampling_threads, size_t evaluation_threads)
{

}

dynamicISALikelihood::~dynamicISALikelihood()
{

}

bool dynamicISALikelihood::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node)
{
	this->varset = varset;
	transformed_variables.setZero(varset->GetNumVariables());

	try {
		boost::property_tree::ptree model_node = likelihood_node.get_child("model");
		std::string model_fn = model_node.get<std::string>("<xmlattr>.file");

		model = std::make_shared<SignalingModel>();
		if (!model->Load(model_fn)) {
			model.reset();
			return false;
		}

		BOOST_FOREACH(const boost::property_tree::ptree::value_type& it, likelihood_node.get_child("")) {
			if (it.first == "experiment") {
				std::unique_ptr<dynamicISAExperiment> experiment = std::make_unique<dynamicISAExperiment>();
				if (!experiment->Load(varset, model, it.second)) {
					return false;
				}
				experiments.push_back(std::move(experiment));
			}
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return true;
}

bool dynamicISALikelihood::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	for (size_t i = 0; i < varset->GetNumVariables(); i++) {
		transformed_variables(i) = varset->TransformVariable(i, values[i]);
	}

	logp = 0.0;
	for (size_t i = 0; i < experiments.size(); i++) {
		Real experiment_logp;
		if (!experiments[i]->EvaluateLogProbability(threadix, transformed_variables, experiment_logp)) {
			return false;
		}
		logp += experiment_logp;
	}

	return true;
}

size_t dynamicISALikelihood::GetNumSignalingMolecules() const
{
	return model->GetNumMolecules();
}

const dynamicISAExperiment* dynamicISALikelihood::GetExperiment(const std::string& experiment) const
{
	for (size_t i = 0; i < experiments.size(); i++) {
		if (experiments[i]->GetName() == experiment) {
			return experiments[i].get();
		}
	}
	LOGERROR("Experiment \"%s\" doesn't exist", experiment.c_str());
	return nullptr;
}
