#include "Utils.h"
#include "dynamicISAExperiment.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"
#include "SignalingModel.h"

dynamicISAExperiment::ModelConditions::ModelConditions()
	: model_ix(std::numeric_limits<size_t>::max())
	, parameter_ix(std::numeric_limits<size_t>::max())
{
}

dynamicISAExperiment::dynamicISAExperiment()
{
}

dynamicISAExperiment::~dynamicISAExperiment()
{
}

bool dynamicISAExperiment::Load(std::shared_ptr<const bcm3::VariableSet> varset, std::shared_ptr<SignalingModel> model, boost::property_tree::ptree xml_node)
{
	bool result = true;

	this->model = model;

	try {
		name = xml_node.get<std::string>("<xmlattr>.name");
		std::string data_fn = xml_node.get<std::string>("<xmlattr>.datafile");
		std::string condition_dim = xml_node.get<std::string>("<xmlattr>.conditions_dim");

		bcm3::NetCDFDataFile file;
		if (!file.Open(data_fn, false)) {
			return false;
		}

		size_t num_conditions;
		result &= file.GetDimensionSize(name, condition_dim, &num_conditions);
		condition_names.resize(num_conditions);
		for (size_t i = 0; i < num_conditions; i++) {
			result &= file.GetValue(name, condition_dim, i, condition_names[i]);
		}

		BOOST_FOREACH(const boost::property_tree::ptree::value_type & var, xml_node.get_child("")) {
			if (var.first == "time_average" || var.first == "timepoint") {
				//if (!ParseDataNode(var.second, other_experiments, &file)) {
				//	return false;
				//}
				Data d;
				d.model_name = var.second.get<std::string>("<xmlattr>.species_name");
				d.model_ix = model->GetMoleculeIxByName(d.model_name);
				if (d.model_ix == std::numeric_limits<size_t>::max()) {
					LOGERROR("Could not find signaling node \"%s\" for data node", d.model_name.c_str());
					return false;
				}

				if (var.first == "timepoint") {
					d.timepoint = var.second.get<Real>("<xmlattr>.time");

					for (size_t i = 0; i < data_timepoints.size(); i++) {
						if (data_timepoints(i) > d.timepoint) {
							LOGERROR("Data must be sorted by increasing time");
							return false;
						}
					}

					data_timepoints.conservativeResize(data_timepoints.size() + 1);
					data_timepoints(data_timepoints.size()-1) = d.timepoint;
				} else {
					d.timepoint = std::numeric_limits<Real>::quiet_NaN();
				}

				std::string data_name = var.second.get<std::string>("<xmlattr>.data_name");
				std::string measurement_dim = var.second.get<std::string>("<xmlattr>.measurement_dim");
				std::string measurement_name = var.second.get<std::string>("<xmlattr>.measurement_name");
				std::string replicate_dim = var.second.get<std::string>("<xmlattr>.replicate_dim", "");

				std::string dimname1, dimname2;
				result &= file.GetDimensionName(name, data_name, 0, dimname1);
				result &= file.GetDimensionName(name, data_name, 1, dimname2);
				if (dimname1 != measurement_dim || dimname2 != condition_dim) {
					LOGERROR("Condition data for \"%s\"/\"%s\" does not have the right dimensions.", measurement_dim.c_str(), measurement_name.c_str());
					return false;
				}
				size_t measurement_dimix = 0;
				result &= file.GetDimensionIx(name, measurement_dim, measurement_name, &measurement_dimix);

				// Check that the variable has the right number of dimensions
				size_t dimcount = 0;
				result &= file.GetDimensionCount(name, data_name, &dimcount);
				if (replicate_dim.empty()) {
					if (dimcount != 2) {
						LOGERROR("Data variable for \"%s\"/\"%s\" has no replicates but is not two-dimensional.", measurement_dim.c_str(), measurement_name.c_str());
						return false;
					}

					d.values.resize(1);
					result &= file.GetValuesDim2(name, data_name, measurement_dimix, 0, condition_names.size(), d.values[0]);
				} else {
					if (dimcount != 3) {
						LOGERROR("Data variable for \"%s\"/\"%s\" has replicates but is not three-dimensional.", measurement_dim.c_str(), measurement_name.c_str());
						return false;
					}
					std::string dimname3;
					result &= file.GetDimensionName(name, data_name, 2, dimname3);
					if (dimname3 != replicate_dim) {
						LOGERROR("Condition data for \"%s\"/\"%s\" does not have the right dimensions.", measurement_dim.c_str(), measurement_name.c_str());
						return false;
					}

					size_t num_replicates;
					result &= file.GetDimensionSize(name, replicate_dim, &num_replicates);
					d.values.resize(num_replicates);
					for (size_t i = 0; i < num_replicates; i++) {
						result &= file.GetValuesDim2(name, data_name, measurement_dimix, 0, i, condition_names.size(), d.values[i]);
					}
				}

				std::string likelihood_param = var.second.get<std::string>("<xmlattr>.likelihood_param");
				d.parameter_base_ix = varset->GetVariableIndex(std::string("base_") + likelihood_param);
				d.parameter_scale_ix = varset->GetVariableIndex(std::string("scale_") + likelihood_param);
				d.parameter_sd_ix = varset->GetVariableIndex(std::string("sd_") + likelihood_param);
				d.parameter_sd_incr_ix = varset->GetVariableIndex(std::string("sdincr_") + likelihood_param);

				data.push_back(d);
			} else if (var.first == "condition") {
				ModelConditions c;
				c.model_name = var.second.get<std::string>("<xmlattr>.species_name");
				c.model_ix = model->GetMoleculeIxByName(c.model_name);
				if (c.model_ix == std::numeric_limits<size_t>::max()) {
					LOGERROR("Could not find signaling node \"%s\" for setting model condition", c.model_name.c_str());
					return false;
				}

				std::string value_str = var.second.get<std::string>("<xmlattr>.value", "");
				std::string data_name = var.second.get<std::string>("<xmlattr>.data_name", "");
				std::string variable_name = var.second.get<std::string>("<xmlattr>.variable_name", "");

				if (!data_name.empty()) {
					// Reference to a variable in the data file; the measurement_dim and name in that dim should also be specified.
					std::string measurement_dim = var.second.get<std::string>("<xmlattr>.measurement_dim");
					std::string measurement_name = var.second.get<std::string>("<xmlattr>.measurement_name");

					// Check that the variable has 2 dimensions
					size_t dimcount = 0;
					result &= file.GetDimensionCount(name, data_name, &dimcount);
					if (dimcount != 2) {
						LOGERROR("Condition data for \"%s\"/\"%s\" is not two-dimensional.", measurement_dim.c_str(), measurement_name.c_str());
						return false;
					}
					std::string dimname1, dimname2;
					result &= file.GetDimensionName(name, data_name, 0, dimname1);
					result &= file.GetDimensionName(name, data_name, 1, dimname2);
					if (dimname1 != measurement_dim || dimname2 != condition_dim) {
						LOGERROR("Condition data for \"%s\"/\"%s\" does not have the right dimensions.", measurement_dim.c_str(), measurement_name.c_str());
						return false;
					}

					size_t measurement_dimix = 0;
					result &= file.GetDimensionIx(name, measurement_dim, measurement_name, &measurement_dimix);
					result &= file.GetValuesDim2(name, data_name, measurement_dimix, 0, condition_names.size(), c.values);
				} else if (!variable_name.empty()) {
					// Condition references an inference variable.
					c.parameter_ix = varset->GetVariableIndex(variable_name, true);
					if (c.parameter_ix == std::numeric_limits<size_t>::max()) {
						LOGERROR("Model condition for variable \"%s\" has specifies a variable reference \"%s\", but can't find it in the variable set", c.model_name.c_str(), variable_name.c_str());
						return false;
					}
				} else if (!value_str.empty()) {
					// Condition is specified directly in the likelihood file
					Real value = boost::lexical_cast<Real>(value_str);
					c.values = VectorReal::Constant(condition_names.size(), value);
				} else {
					LOGERROR("Model condition for variable \"%s\" has neither data_name, variable_name or value specified", c.model_name.c_str());
					return false;
				}

				condition_specifications.push_back(c);
			} else if (var.first == "expression_level") {
#if 0
				ExpressionLevel el;
				el.model_name = var.second.get<std::string>("<xmlattr>.species_name");
				el.model_ix = network->GetSignalingMoleculeIxByName(el.model_name);
				if (el.model_ix == std::numeric_limits<size_t>::max()) {
					LOGERROR("Could not find signaling node \"%s\" for setting expression level", el.model_name.c_str());
					return false;
				}

				std::string value_str = var.second.get<std::string>("<xmlattr>.value", "");
				std::string data_name = var.second.get<std::string>("<xmlattr>.data_name", "");

				if (!data_name.empty()) {
					std::string data_name_without_index;
					std::vector<size_t> data_indices;
					if (!ParseDataFileReference(data_name, data_name_without_index, data_indices, &file)) {
						LOGERROR("Unable to parse data reference \"%s\" for expression level for species \"%s\"", data_name.c_str(), el.model_name.c_str());
						return false;
					}
					if (data_indices.size() != 2) {
						LOGERROR("Expression level for species \"%s\" specifies a data reference \"%s\", but it does not have 2 dimensions (should be: [genes/epitopes,cell_lines])", el.model_name.c_str(), data_name.c_str());
						return false;
					}

					// TODO - verify dimension 2 is cell lines

					result &= file.GetValuesDim2(Name, data_name_without_index, data_indices[0], 0, cell_lines.size(), el.values);
				} else if (!value_str.empty()) {
					Real value = boost::lexical_cast<Real>(value_str);
					el.values = VectorReal::Constant(cell_lines.size(), value);
				} else {
					LOGERROR("Expression level for variable \"%s\" has neither data_name nor value specified", el.model_name.c_str());
					return false;
				}

				std::string base_parameter_str = var.second.get<std::string>("<xmlattr>.base_parameter", std::string("base_expression[") + el.model_name + std::string("]"));
				el.base_parameter_ix = varset->GetVariableIndex(base_parameter_str, false);
				std::string scale_parameter_str = var.second.get<std::string>("<xmlattr>.scale_parameter", std::string("scale_expression[") + el.model_name + std::string("]"));
				el.scale_parameter_ix = varset->GetVariableIndex(scale_parameter_str, false);

				expression_levels.push_back(el);
#endif
			} else if (var.first == "treatment") {
				std::string species_name = var.second.get<std::string>("<xmlattr>.species_name");
				size_t model_ix = model->GetMoleculeIxByName(species_name);
				if (model_ix == std::numeric_limits<size_t>::max()) {
					LOGERROR("Could not find signaling node \"%s\" for treatment node", species_name.c_str());
					return false;
				}

				Real time = var.second.get<Real>("<xmlattr>.time");
				Real value = var.second.get<Real>("<xmlattr>.value");

				treatments.push_back(std::tuple<size_t, Real, Real>(model_ix, time, value));
			}
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	base_param_values.setConstant(model->GetNumMolecules(), 0.0);
	decay_param_values.setConstant(model->GetNumMolecules(), 0.0);
	strength_param_values.setConstant(model->GetNumMolecules(), 0.0);
	inhib_param_values.setConstant(model->GetNumMolecules(), 0.0);
	expression_mixing_param_values.setConstant(model->GetNumMolecules(), 1.0);
	fixed_species_values.setConstant(model->GetNumMolecules(), 0.0);
	expression_levels.setConstant(model->GetNumMolecules(), 1.0);

	// Map inference variables
	base_param_map.resize(model->GetNumMolecules());
	for (size_t i = 0; i < model->GetNumMolecules(); i++) {
		base_param_map[i] = varset->GetVariableIndex(std::string("base_") + model->GetMoleculeName(i), false);
		if (base_param_map[i] == std::numeric_limits<size_t>::max()) {
			// Has this been specified by a condition?
			bool specified = false;
			for (auto csi = condition_specifications.begin(); csi != condition_specifications.end(); ++csi) {
				if (csi->model_ix == i) {
					specified = true;
				}
			}
			if (!specified) {
				LOGWARNING("No base value specification for signaling molecule \"%s\" in experiment \"%s\", assuming 0", model->GetMoleculeName(i).c_str(), name.c_str());
			}
		}
	}
	decay_param_map.resize(model->GetNumMolecules());
	for (size_t i = 0; i < model->GetNumMolecules(); i++) {
		decay_param_map[i] = varset->GetVariableIndex(std::string("decay_") + model->GetMoleculeName(i), false);
		if (decay_param_map[i] == std::numeric_limits<size_t>::max()) {
			// Has this been specified by a condition?
			bool specified = false;
			for (auto csi = condition_specifications.begin(); csi != condition_specifications.end(); ++csi) {
				if (csi->model_ix == i) {
					specified = true;
				}
			}
			if (!specified) {
				LOGWARNING("No decay value specification for signaling molecule \"%s\" in experiment \"%s\", assuming 0", model->GetMoleculeName(i).c_str(), name.c_str());
			}
		}
	}
	strength_param_map.resize(model->GetNumReactions());
	for (size_t i = 0; i < model->GetNumReactions(); i++) {
		strength_param_map[i] = varset->GetVariableIndex(model->GetStrengthName(i));
		if (strength_param_map[i] == std::numeric_limits<size_t>::max()) {
			LOGERROR("No strength value specification for signaling link \"%s\" in experiment \"%s\"", model->GetStrengthName(i).c_str(), name.c_str());
			return false;
		}
	}

	modeled_activities.setConstant(condition_names.size(), model->GetNumMolecules(), std::numeric_limits<Real>::quiet_NaN());
	modeled_data.resize(data.size(), VectorReal::Constant(condition_names.size(), std::numeric_limits<Real>::quiet_NaN()));

	return result;
}

bool dynamicISAExperiment::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	for (size_t i = 0; i < condition_names.size(); i++) {
		Real condition_logp = 0.0;

		for (size_t i = 0; i < model->GetNumMolecules(); i++) {
			if (base_param_map[i] != std::numeric_limits<size_t>::max()) {
				base_param_values(i) = values(base_param_map[i]);
			} else {
				base_param_values(i) = 0.0;
			}
			if (decay_param_map[i] != std::numeric_limits<size_t>::max()) {
				decay_param_values(i) = values(decay_param_map[i]);
			} else {
				decay_param_values(i) = 0.0;
			}
		}
		for (size_t i = 0; i < model->GetNumReactions(); i++) {
			strength_param_values[i] = values(strength_param_map[i]);
		}

		fixed_species_values.setConstant(model->GetNumMolecules(), std::numeric_limits<Real>::quiet_NaN());
		for (auto csi = condition_specifications.begin(); csi != condition_specifications.end(); ++csi) {
			if (csi->parameter_ix != std::numeric_limits<size_t>::max()) {
				fixed_species_values(csi->model_ix) = values(csi->parameter_ix);
			} else {
				fixed_species_values(csi->model_ix) = csi->values(i);
			}
		}

#if 1
		VectorReal activities;
		bool result = model->Calculate(base_param_values,
			decay_param_values,
			strength_param_values,
			inhib_param_values,
			expression_mixing_param_values,
			fixed_species_values,
			expression_levels,
			24.0 * 3600.0,
			activities);
		if (!result) {
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}

		modeled_activities.row(i) = activities;

		for (size_t j = 0; j < data.size(); j++) {
			Data& d = data[j];
			for (size_t k = 0; k < d.values.size(); k++) {
				Real x = values(d.parameter_base_ix) + values(d.parameter_scale_ix) * activities(d.model_ix);
				modeled_data[j](i) = x;

				Real y = d.values[k](i);
				Real sd = values(d.parameter_sd_ix) + fabs(x) * values(d.parameter_sd_incr_ix);
				if (!std::isnan(y)) {
					condition_logp += bcm3::LogPdfTnu3(y, x, sd);
				}
			}
		}
#else
		std::vector<VectorReal> activities(data_timepoints.size(), VectorReal::Constant(model->GetNumMolecules(), std::numeric_limits<Real>::quiet_NaN()));
		bool result = model->Calculate(base_param_values,
			decay_param_values,
			strength_param_values,
			inhib_param_values,
			expression_mixing_param_values,
			fixed_species_values,
			expression_levels,
			data_timepoints,
			treatments,
			activities);
		if (!result) {
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}

		for (size_t j = 0; j < data.size(); j++) {
			Data& d = data[j];
			for (size_t k = 0; k < d.values.size(); k++) {
				Real x = values(d.parameter_base_ix) + values(d.parameter_scale_ix) * activities[j](d.model_ix);
				modeled_data[j](i) = x;

				Real y = d.values[k](i);
				Real sd = values(d.parameter_sd_ix) + fabs(x) * values(d.parameter_sd_incr_ix);
				if (!std::isnan(y)) {
					condition_logp += bcm3::LogPdfTnu3(y, x, sd);
				}
			}
		}
#endif

		logp += condition_logp;
	}

	return true;
}

void dynamicISAExperiment::GetObservedData(size_t data_ix, MatrixReal& out_values) const
{
	const Data& d = data[data_ix];
	out_values.resize(condition_names.size(), d.values.size());
	for (size_t k = 0; k < d.values.size(); k++) {
		out_values.col(k) = d.values[k];
	}
}

void dynamicISAExperiment::GetModeledData(size_t data_ix, VectorReal& out_values) const
{
	out_values = modeled_data[data_ix];
}

void dynamicISAExperiment::GetModeledActivities(MatrixReal& out_values) const
{
	out_values = modeled_activities;
}
