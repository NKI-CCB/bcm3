#include "Utils.h"
#include "PriorIndependence.h"
#include "UnivariateMarginal.h"
#include "VariableSet.h"

#include <boost/property_tree/xml_parser.hpp>

namespace bcm3 {

PriorIndependence::PriorIndependence(std::shared_ptr<VariableSet> varset)
	: Prior(varset)
{
	Variables.resize(varset->GetNumVariables());
}

PriorIndependence::~PriorIndependence()
{
}

bool PriorIndependence::LoadFromXML(const boost::property_tree::ptree& xml_node)
{
	try {
		size_t variable_ix = 0;

		BOOST_FOREACH(const boost::property_tree::ptree::value_type& var, xml_node.get_child("")) {
			if (var.first == "variable") {
				std::string name = var.second.get<std::string>("<xmlattr>.name");

				auto m = std::make_unique<UnivariateMarginal>();
				if (!m->Initialize(var.second)) {
					return false;
				}
				if (!SetMarginal(name, m)) {
					return false;
				}

				variable_ix++;
			}
		}

		return true;
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing prior file: %s", e.what());
		return false;
	}
}

bool PriorIndependence::SetMarginal(const std::string& name, std::unique_ptr<UnivariateMarginal>& marginal)
{
	size_t ix = varset->GetVariableIndex(name);
	if (ix == std::numeric_limits<size_t>::max()) {
		LOGERROR("Attempt to set marginal for unknown variable \"%s\"", name.c_str());
		return false;
	}

	if (Variables[ix]) {
		LOGERROR("Attempt to set marginal for variable \"%s\" twice", name.c_str());
		return false;
	}

	LowerBounds(ix) = marginal->GetLowerBound();
	UpperBounds(ix) = marginal->GetUpperBound();
	Variables[ix].swap(marginal);
	return true;
}

bool PriorIndependence::EvaluateLogPDF(size_t threadix, const VectorReal& values, Real& logp) const
{
	logp = 0.0;

	for (size_t i = 0; i < Variables.size(); i++) {
		ASSERT(Variables[i]);

		Real varlogp;
		if (!Variables[i]->EvaluateLogPDF(varlogp, values(i))) {
			LOG("Error evaluating prior for variable %d", i);
			return false;
		}
		if (varlogp != varlogp) {
			LOG("NaN in evaluating prior for variable %d, value=%g", i, values[i]);
		}
		logp += varlogp;
	}

	return true;
}

bool PriorIndependence::Sample(VectorReal& values, RNG* rng) const
{
	values.resize(Variables.size());

	for (size_t i = 0; i < Variables.size(); i++) {
		ASSERT(Variables[i]);

		if (!Variables[i]->Sample(values(i), rng)) {
			return false;
		}
	}

	return true;
}

bool PriorIndependence::EvaluateMarginalMean(size_t i, Real& mean) const
{
	if (i >= varset->GetNumVariables()) {
		return false;
	}
	ASSERT(Variables[i]);
	return Variables[i]->EvaluateMean(mean);
}

bool PriorIndependence::EvaluateMarginalVariance(size_t i, Real& var) const
{
	if (i >= varset->GetNumVariables()) {
		return false;
	}
	ASSERT(Variables[i]);
	return Variables[i]->EvaluateVariance(var);
}

}
