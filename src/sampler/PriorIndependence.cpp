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

					bool multivariate = var.second.get<bool>("<xmlattr>.multivariate", false);

					if (multivariate) {
						size_t id = var.second.get<size_t>("<xmlattr>.id");
						if (id == 0) {
							LOGERROR("Multivariate distribution IDs should start at 1.");
							return false;
						} else if (id <= MultivariateDistributions.size()) {
						} else {
							MultivariateDistributions.resize(id);
							MultivariateDistributions[id - 1] = std::make_unique<MultivariateMarginal>();
							// Only dirichlet supported for now.
							if (var.second.get<std::string>("<xmlattr>.distribution") != "dirichlet") {
								LOGERROR("Multivariate distribution of unknown type (only dirichlet supported).");
								return false;
							}
							MultivariateDistributions[id - 1]->CreateDirichlet();
							MultivariateVarIx.resize(id);
							MultivariateVarIx[id - 1] = variable_ix;
						}

						Real alpha = var.second.get<Real>("<xmlattr>.alpha");
						size_t mv_var_ix;
						if (MultivariateDistributions[id - 1]->AddVariable(alpha, mv_var_ix)) {
							MultivariateMembership.push_back(id - 1);
							if (variable_ix != MultivariateVarIx[id - 1] + MultivariateDistributions[id - 1]->GetSize() - 1) {
								LOGERROR("All variables in a multivariate distribution should follow each other directly");
								return false;
							}
						} else {
							return false;
						}

						LowerBounds[variable_ix] = MultivariateDistributions[id - 1]->GetLowerBound(mv_var_ix);
						UpperBounds[variable_ix] = MultivariateDistributions[id - 1]->GetUpperBound(mv_var_ix);

						size_t ix = varset->GetVariableIndex(name);
						Variables[ix].reset();
					} else {
						auto m = std::make_unique<UnivariateMarginal>();
						if (!m->Initialize(var.second)) {
							return false;
						}
						if (!SetMarginal(name, m)) {
							return false;
						}

						MultivariateMembership.push_back(std::numeric_limits<size_t>::max());
					}

					variable_ix++;
				}
			}

			for (size_t i = 0; i < MultivariateDistributions.size(); i++) {
				if (!MultivariateDistributions[i]->Initialize()) {
					return false;
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

		for (size_t i = 0; i < MultivariateDistributions.size(); i++) {
			Real mvlogp;
			MultivariateDistributions[i]->EvaluateLogPDF(mvlogp, &values[MultivariateVarIx[i]]);
			logp += mvlogp;
		}

		for (size_t i = 0; i < Variables.size(); i++) {
			if (Variables[i] == NULL) {
				// Must be part of an mvdist
				continue;
			}

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

		for (size_t i = 0; i < MultivariateDistributions.size(); i++) {
			MultivariateDistributions[i]->Sample(&values[MultivariateVarIx[i]], rng);
		}

		for (size_t i = 0; i < Variables.size(); i++) {
			if (Variables[i] == NULL) {
				// Must be part of an mvdist
				continue;
			}

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

		if (!Variables[i]) {
			size_t mvdistix = MultivariateMembership.at(i);
			return MultivariateDistributions[mvdistix]->EvaluateMarginalMean(i - MultivariateVarIx[mvdistix], mean);
		} else {
			return Variables[i]->EvaluateMean(mean);
		}
	}

	bool PriorIndependence::EvaluateMarginalVariance(size_t i, Real& var) const
	{
		if (i >= varset->GetNumVariables()) {
			return false;
		}

		if (!Variables[i]) {
			size_t mvdistix = MultivariateMembership.at(i);
			return MultivariateDistributions[mvdistix]->EvaluateMarginalVariance(i - MultivariateVarIx[mvdistix], var);
		} else {
			return Variables[i]->EvaluateVariance(var);
		}
	}

	size_t PriorIndependence::GetNumDirichletDistributions() const
	{
		// All multivariate distributions are dirichlets for now
		return MultivariateDistributions.size();
	}

	void PriorIndependence::GetDirichletMembership(ptrdiff_t ix, ptrdiff_t& residual_ix, std::vector<ptrdiff_t>& var_ixs) const
	{
		var_ixs.clear();
		for (ptrdiff_t i = 0; i < Variables.size(); i++) {
			if (MultivariateMembership.at(i) == ix) {
				var_ixs.push_back(i);
			}
		}
		residual_ix = *var_ixs.rbegin();
	}

}
