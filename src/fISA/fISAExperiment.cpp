#include "Utils.h"
#include "fISAExperiment.h"
#include "fISAExperimentDrugRange.h"
#include "fISAExperimentIncucyteSequential.h"
#include "fISAExperimentSingleCondition.h"
#include "SignalingNetwork.h"
#include "ProbabilityDistributions.h"

fISAExperiment::ModelConditions::ModelConditions()
	: model_ix(std::numeric_limits<size_t>::max())
	, parameter_ix(std::numeric_limits<size_t>::max())
{
}

fISAExperiment::ExpressionLevel::ExpressionLevel()
	: model_ix(std::numeric_limits<size_t>::max())
	, base_parameter_ix(std::numeric_limits<size_t>::max())
	, scale_parameter_ix(std::numeric_limits<size_t>::max())
{
}

fISAExperiment::DataPartBase::DataPartBase()
	: likelihood_fn(LF_Invalid)
	, use_base(false)
	, use_scale(false)
	, fixed_base(false)
	, fixed_sd(false)
	, weight(1.0)
	, model_ix(std::numeric_limits<size_t>::max())
	, parameter_base_ix(std::numeric_limits<size_t>::max())
	, parameter_scale_ix(std::numeric_limits<size_t>::max())
	, parameter_sd_ix(std::numeric_limits<size_t>::max())
	, fixed_sd_value(0.0)
	, fixed_sd_scale_value(0.0)
	, data_is_inactive_form(false)
	, expression_ix(std::numeric_limits<size_t>::max())
	, expression_ix2(std::numeric_limits<size_t>::max())
	, expression_sum_parameter_ix(std::numeric_limits<size_t>::max())
	, model_ix2(std::numeric_limits<size_t>::max())
	, data_sum_parameter_ix(std::numeric_limits<size_t>::max())
	, parameter_sd_base_ix(std::numeric_limits<size_t>::max())
{
}

fISAExperiment::fISAExperiment()
	: network(nullptr)
	, varset(nullptr)
	, task_manager(nullptr)
	//, safe_bayes_ix(std::numeric_limits<size_t>::max())
{
}

fISAExperiment::~fISAExperiment()
{
}

std::unique_ptr<fISAExperiment> fISAExperiment::Create(const boost::property_tree::ptree& xml_node, SignalingNetwork* network, std::shared_ptr<const bcm3::VariableSet> varset, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, size_t evaluation_threads, bcm3::TaskManager* task_manager)
{
	std::unique_ptr<fISAExperiment> experiment;

	try {
		std::string type = xml_node.get<std::string>("<xmlattr>.type", "single_condition");

		if (type == "single_condition") {
			experiment = std::make_unique<fISAExperimentSingleCondition>();
		}
#if TODO
		else if (type == "drug_range") {
			experiment = std::make_unique<fISAExperimentDrugRange>(numthreads);
		} else if (type == "incucyte_sequential") {
			experiment = std::make_unique<fISAExperimentIncucyteSequential>(numthreads);
		}
#endif
		else {
			LOGERROR("Unknown experiment type \"%s\"", type.c_str());
			return NULL;
		}

		experiment->Name = xml_node.get<std::string>("<xmlattr>.name");
		experiment->DataFile = xml_node.get<std::string>("<xmlattr>.datafile");
		experiment->network = network;
		experiment->varset = varset;
		experiment->task_manager = task_manager;

		if (!experiment->Load(xml_node, other_experiments)) {
			experiment.reset();
		}

	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		experiment.reset();
	}

	if (experiment) {
		experiment->transformed_values.resize(varset->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());
		if (!experiment->InitializeParallelData(evaluation_threads)) {
			experiment.reset();
		}
	}

	return experiment;
}

#if 0
void fISAExperiment::SetBootstrap(unsigned long rngseed)
{
	bcm::RNG_SFMT rng(rngseed);

	cell_lines_bootstrap_count.resize(cell_lines.size(), 0);
	for (size_t i = 0; i < cell_lines.size(); i++) {
		unsigned long s = rng.GetUniformInt(cell_lines.size());
		cell_lines_bootstrap_count[s]++;
	}
}

void fISAExperiment::SetSafeBayesIx(size_t ix)
{
	safe_bayes_ix = ix;
}

void fISAExperiment::AddSafeBayesCvIx(size_t ix)
{
	safe_bayes_cv_ix.insert(ix);
}
#endif

bool fISAExperiment::Load(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	bool result = true;
	bcm3::NetCDFDataFile file;
	
	if (!file.Open(DataFile, false)) {
		return false;
	}

	size_t num_cell_lines;
	result &= file.GetDimensionSize(Name, "cell_lines", &num_cell_lines);
	cell_lines.resize(num_cell_lines);
	for (size_t i = 0; i < cell_lines.size(); i++) {
		result &= file.GetValue(Name, "cell_lines", i, cell_lines[i]);
	}

	if (!LoadTypeSpecificNodes(xml_node, file)) {
		return false;
	}
	
	try {
		BOOST_FOREACH(const boost::property_tree::ptree::value_type& var, xml_node.get_child("")) {
			if (var.first == "data") {
				if (!ParseDataNode(var.second, other_experiments, file)) {
					return false;
				}
			} else if (var.first == "mutation" || var.first == "condition") {
				ModelConditions c;
				c.model_name = var.second.get<std::string>("<xmlattr>.species_name");
				c.model_ix = network->GetSignalingMoleculeIxByName(c.model_name);
				if (c.model_ix == std::numeric_limits<size_t>::max()) {
					LOGERROR("Could not find signaling node \"%s\" for setting model condition", c.model_name.c_str());
					return false;
				}

				std::string value_str = var.second.get<std::string>("<xmlattr>.value", "");
				std::string data_name = var.second.get<std::string>("<xmlattr>.data_name", "");
				std::string variable_name = var.second.get<std::string>("<xmlattr>.variable_name", "");

				if (!data_name.empty()) {
					std::string data_name_without_index;
					std::vector<size_t> data_indices;
					if (!ParseDataFileReference(data_name, data_name_without_index, data_indices, file)) {
						LOGERROR("Unable to parse data reference \"%s\" for model condition for species \"%s\" in experiment \"%s\"", data_name.c_str(), c.model_name.c_str(), Name.c_str());
						return false;
					}
					if (data_indices.size() != 2) {
						LOGERROR("Model condition for species \"%s\" specifies a data reference \"%s\", but it does not have 2 dimensions (should be: [genes/epitopes,cell_lines])", c.model_name.c_str(), data_name.c_str());
						return false;
					}

					// TODO - verify dimension 2 is cell lines

					result &= file.GetValuesDim2(Name, data_name_without_index, data_indices[0], 0, cell_lines.size(), c.values);
				} else if (!variable_name.empty()) {
					c.parameter_ix = varset->GetVariableIndex(variable_name, true);
					if (c.parameter_ix == std::numeric_limits<size_t>::max()) {
						LOGERROR("Model condition for variable \"%s\" has specifies a variable reference \"%s\", but can't find it in the variable set", c.model_name.c_str(), variable_name.c_str());
						return false;
					}
				} else if (!value_str.empty()) {
					Real value = boost::lexical_cast<Real>(value_str);
					c.values = VectorReal::Constant(cell_lines.size(), value);
				} else {
					LOGERROR("Model condition for variable \"%s\" has neither data_name, variable_name or value specified", c.model_name.c_str());
					return false;
				}

				conditions.push_back(c);
			} else if (var.first == "expression_level") {
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
					if (!ParseDataFileReference(data_name, data_name_without_index, data_indices, file)) {
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
			}
		}
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	file.Close();
	return true;
}

bool fISAExperiment::ParseDataPartBase(const boost::property_tree::ptree& xml_node, DataPartBase& p, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	p.model_name = xml_node.get<std::string>("<xmlattr>.species_name");
	p.data_name  = xml_node.get<std::string>("<xmlattr>.data_name");
	p.base_scale_sd_suffix  = xml_node.get<std::string>("<xmlattr>.base_scale_sd_suffix", "");

	std::string likelihood_fn_str = xml_node.get<std::string>("<xmlattr>.likelihood_function", "studentt");
	if (likelihood_fn_str == "normal") {
		p.likelihood_fn = LF_Normal;
	} else if (likelihood_fn_str == "truncated_normal") {
		p.likelihood_fn = LF_TruncatedNormal;
	}  else if (likelihood_fn_str == "studentt") {
		p.likelihood_fn = LF_StudentT;
	} else if (likelihood_fn_str == "truncated_t") {
		p.likelihood_fn = LF_TruncatedStudentT;
	} else if (likelihood_fn_str == "binomial") {
		p.likelihood_fn = LF_Binomial;
	} else if (likelihood_fn_str == "beta") {
		p.likelihood_fn = LF_Beta;
	} else if (likelihood_fn_str == "ordered_probit") {
		p.likelihood_fn = LF_OrderedProbit;
	} else if (likelihood_fn_str == "ordered_robit") {
		p.likelihood_fn = LF_OrderedRobit;
	} else {
		LOGERROR("Unknown likelihood function \"%s\"", likelihood_fn_str.c_str());
		return false;
	}

	p.use_base = xml_node.get<bool>("<xmlattr>.use_base", true);
	p.use_scale = xml_node.get<bool>("<xmlattr>.use_scale", true);
	p.scale_var_with_mean = xml_node.get<bool>("<xmlattr>.scale_var_with_mean", true);
	p.fixed_sd_scale_value = xml_node.get<Real>("<xmlattr>.sd_scale", std::numeric_limits<Real>::quiet_NaN());
	p.weight = xml_node.get<Real>("<xmlattr>.weight", 1.0);

	p.data_is_inactive_form = xml_node.get<bool>("<xmlattr>.data_is_inactive_form", false);
	p.scale_per_cell_line = xml_node.get<bool>("<xmlattr>.scale_per_cell_line", false);

	if (p.use_base) {
		std::string base_str = xml_node.get<std::string>("<xmlattr>.base", "");
		if (base_str.empty()) {
			base_str = std::string("base_") + p.base_scale_sd_suffix;
		}

		// Is the base a parameter in the prior?
		p.parameter_base_ix = varset->GetVariableIndex(base_str, false);
		if (p.parameter_base_ix == std::numeric_limits<size_t>::max()) {
			// No, it isn't - try to cast it to a Real
			try {
				p.fixed_base_value = boost::lexical_cast<Real>(base_str);
				p.fixed_base = true;
			}
			catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for data base parameter \"%s\", and could also not cast it to a constant real value: %s", base_str.c_str(), e.what());
				return false;
			}
		}
		else {
			p.fixed_base = false;
		}
	}
				
	if (p.use_scale) {
		p.parameter_scale_ix = varset->GetVariableIndex(std::string("scale_") + p.base_scale_sd_suffix);
		if (p.parameter_scale_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Could not find variable for data scale parameter with suffix \"%s\"", p.base_scale_sd_suffix.c_str());
			return false;
		}
	}

	if (p.likelihood_fn == LF_Normal || p.likelihood_fn == LF_TruncatedNormal || p.likelihood_fn == LF_StudentT || p.likelihood_fn == LF_TruncatedStudentT) {
		std::string sd_str = xml_node.get<std::string>("<xmlattr>.sd", "");
		if (sd_str.empty()) {
			sd_str = std::string("sd_") + p.base_scale_sd_suffix;
		}

		// Is the standard deviation a parameter in the prior?
		p.parameter_sd_ix = varset->GetVariableIndex(sd_str, false);
		if (p.parameter_sd_ix == std::numeric_limits<size_t>::max()) {
			// No, it isn't - try to cast it to a Real
			try {
				p.fixed_sd_value = boost::lexical_cast<Real>(sd_str);
				p.fixed_sd = true;
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for data standard deviation parameter \"%s\", and could also not cast it to a constant real value: %s", sd_str.c_str(), e.what());
				return false;
			}
		}
		else {
			p.fixed_sd = false;
		}

		p.parameter_sd_base_ix = varset->GetVariableIndex(sd_str + "_base", false);
	}

	if (p.likelihood_fn == LF_Beta) {
		std::string precision_str = xml_node.get<std::string>("<xmlattr>.precision", "");
		if (precision_str.empty()) {
			precision_str = std::string("precision_") + p.base_scale_sd_suffix;
		}

		// Is the variance a parameter in the prior?
		p.parameter_precision_ix = varset->GetVariableIndex(precision_str, false);
		if (p.parameter_precision_ix == std::numeric_limits<size_t>::max()) {
			// No, it isn't - try to cast it to a Real
			try {
				p.fixed_precision_value = boost::lexical_cast<Real>(precision_str);
				p.fixed_sd = true;
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for data precision parameter \"%s\", and could also not cast it to a constant real value: %s", precision_str.c_str(), e.what());
				return false;
			}
		} else {
			p.fixed_sd = false;
		}
	}

	if (p.likelihood_fn == LF_OrderedProbit || p.likelihood_fn == LF_OrderedRobit) {
		size_t categories = xml_node.get<size_t>("<xmlattr>.categories");
		p.ordered_probit_thresholds.resize(categories-1, 0);

#if 1
		std::string threshold_str = xml_node.get<std::string>("<xmlattr>.thresholds", "");
		if (!threshold_str.empty()) {
			for (size_t i = 0; i < categories - 1; i++) {
				std::string varstr = threshold_str + "[" + std::to_string(i) + "]";
				size_t var_ix = varset->GetVariableIndex(varstr, false);
				if (p.parameter_precision_ix == std::numeric_limits<size_t>::max()) {
					LOGERROR("Could not find variable for ordered probit threshold \"%s\"", varstr.c_str());
					return false;
				}
				p.ordered_probit_thresholds[i] = var_ix;
			}
		}
#endif
	}

	// Is the data a sum of two model variables?
	if (p.model_name.find('+') != std::string::npos) {
		std::vector<std::string> vars;
		bcm3::tokenize(p.model_name, vars, "+");
		if (vars.size() != 2) {
			LOGERROR("Species name in data reference \"%s\" has a '+' sign, indicating the data is a sum of model species; but there are %z species and only 2 species are supported.", p.model_name.c_str(), vars.size());
			return false;
		}
		p.model_ix = network->GetSignalingMoleculeIxByName(vars[0]);
		if (p.model_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Could not find signaling node \"%s\"", vars[0].c_str());
			return false;
		}

		p.model_ix2 = network->GetSignalingMoleculeIxByName(vars[1]);
		if (p.model_ix2 == std::numeric_limits<size_t>::max()) {
			LOGERROR("Could not find signaling node \"%s\"", vars[1].c_str());
			return false;
		}

		std::string data_sum_parameter_str = xml_node.get<std::string>("<xmlattr>.data_sum_parameter", "");
		if (!data_sum_parameter_str.empty()) {
			p.data_sum_parameter_ix = varset->GetVariableIndex(data_sum_parameter_str, true);
		} else {
			p.data_sum_parameter_ix = std::numeric_limits<size_t>::max();
		}
	} else {
		p.model_ix = network->GetSignalingMoleculeIxByName(p.model_name);
		if (p.model_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Could not find signaling node \"%s\"", p.model_name.c_str());
			return false;
		}

		p.model_ix2 = std::numeric_limits<size_t>::max();
		p.data_sum_parameter_ix = std::numeric_limits<size_t>::max();
	}

	std::string type = xml_node.get<std::string>("<xmlattr>.type", "");
	if (type == "relative") {
		p.relative_to_experiment = xml_node.get<std::string>("<xmlattr>.relative_reference");
					
		// Make sure the other experiment is defined before this one
		bool found = false;
		for (std::vector< std::unique_ptr<fISAExperiment> >::const_iterator ei = other_experiments.begin(); ei != other_experiments.end(); ++ei) {
			if ((*ei)->Name == p.relative_to_experiment) {
				found = true;
				break;
			}
		}
		if (!found) {
			LOGERROR("A data point for node \"%s\" is specified as relative to experiment \"%s\", but that experiment has not been defined. Make sure relative experiments are defined in the right order.", p.data_name.c_str(), p.relative_to_experiment.c_str());
			return false;
		}
	}

	std::string expression_str = xml_node.get<std::string>("<xmlattr>.expression", "");
	if (!expression_str.empty()) {
		// Is the expression a sum of two model variables?
		if (expression_str.find('+') != std::string::npos) {
			std::vector<std::string> vars;
			bcm3::tokenize(expression_str, vars, "+");
			if (vars.size() != 2) {
				LOGERROR("Expression in data reference \"%s\" has a '+' sign, indicating the data is a sum of model species; but there are %z species and only 2 species are supported.", p.model_name.c_str(), vars.size());
				return false;
			}
			p.expression_ix = network->GetSignalingMoleculeIxByName(vars[0]);
			if (p.expression_ix == std::numeric_limits<size_t>::max()) {
				LOGERROR("Expression has been specified for node \"%s\" as \"%s\", but that molecule is not found \"%s\" in the signaling network", p.data_name.c_str(), expression_str.c_str(), vars[0].c_str());
				return false;
			}

			p.expression_ix2 = network->GetSignalingMoleculeIxByName(vars[1]);
			if (p.expression_ix2 == std::numeric_limits<size_t>::max()) {
				LOGERROR("Expression has been specified for node \"%s\" as \"%s\", but that molecule is not found \"%s\" in the signaling network", p.data_name.c_str(), expression_str.c_str(), vars[1].c_str());
				return false;
			}

			std::string expression_sum_parameter_str = xml_node.get<std::string>("<xmlattr>.expression_sum_parameter", "");
			if (!expression_sum_parameter_str.empty()) {
				p.expression_sum_parameter_ix = varset->GetVariableIndex(expression_sum_parameter_str, true);
			} else {
				p.expression_sum_parameter_ix = std::numeric_limits<size_t>::max();
			}
		} else {
			p.expression_ix = network->GetSignalingMoleculeIxByName(expression_str);
			if (p.expression_ix == std::numeric_limits<size_t>::max()) {
				LOGERROR("Expression has been specified for node \"%s\" as \"%s\", but that molecule is not found in the signaling network", p.data_name.c_str(), expression_str.c_str());
				return false;
			}
		}
	}

	return true;
}

void fISAExperiment::PrepareActivitiesCalculation(VectorReal& activities, VectorReal& expression, const Real* transformed_values, size_t cell_ix) const
{
	// Initialize activities to NaN, so that SignalingNetwork knows which values it still has to calculate.
	activities.setConstant(std::numeric_limits<Real>::quiet_NaN());

	// Clear all drug concentrations and set the new concentration for this drug
	network->UpdateActivitiesForDrugSpecies(activities, 0.0);

	// Fill the activities for the nodes which are set to a value in the model conditions.
	for (size_t mci = 0; mci < conditions.size(); mci++) {
		const ModelConditions& c = conditions[mci];
		if (c.parameter_ix != std::numeric_limits<size_t>::max()) {
			activities(c.model_ix) = transformed_values[c.parameter_ix];
		} else {
			activities(c.model_ix) = c.values(cell_ix);
		}
	}
	
	// Calculate expression values
	expression = VectorReal::Ones(network->GetMoleculeCount());
	for (std::vector<ExpressionLevel>::const_iterator eli = expression_levels.begin(); eli != expression_levels.end(); ++eli) {
		if (eli->base_parameter_ix != std::numeric_limits<size_t>::max()) {
			if (eli->scale_parameter_ix != std::numeric_limits<size_t>::max()) {
				const Real eb = transformed_values[eli->base_parameter_ix];
				const Real es = transformed_values[eli->scale_parameter_ix];
				Real e = (eli->values(cell_ix) - eb) / es;
				e = std::max(std::min(e, 1.0), 0.0);
				//e = bcm3::logistic(e);
				expression(eli->model_ix) = e;
			} else {
				const Real eb = transformed_values[eli->base_parameter_ix];
				Real e = (eli->values(cell_ix) - eb) / (1.0 - eb);
				e = std::max(std::min(e, 1.0), 0.0);
				//e = bcm3::logistic(e);
				expression(eli->model_ix) = e;
			}
		} else {
			Real e = eli->values(cell_ix);
			//e = bcm3::logistic(e);
			expression(eli->model_ix) = e;
		}
	}
}

bool fISAExperiment::ParseDataFileReference(const std::string& data_name, std::string& data_name_without_index, std::vector<size_t>& data_indices, const bcm3::NetCDFDataFile& datafile)
{
	//
	// The data name parsing is copied from sbmlodeinf/LikelihoodPart.cpp
	//

	data_indices.clear();

	size_t indexing_start = data_name.find_first_of('[');
	if (indexing_start == data_name.npos) {
		// No indexing specified
		data_name_without_index = data_name;
		return true;
	}

	size_t indexing_end = data_name.find_first_of(']');
	if (indexing_end == data_name.npos) {
		LOGERROR("Syntax error in data name: missing closing bracket for indexing");
		return false;
	}

	data_name_without_index = data_name.substr(0, indexing_start);

	std::string indexingstr = data_name.substr(indexing_start+1, indexing_end - indexing_start - 1);
	std::vector<std::string> indexing_parts;
	bcm3::tokenize(indexingstr, indexing_parts, ",");

	for (size_t i = 0; i < indexing_parts.size(); i++) {
		if (indexing_parts[i] == ":") {
			data_indices.push_back(std::numeric_limits<size_t>::max());
		} else {
			if (indexing_parts[i].front() == '\'' && indexing_parts[i].back() == '\'') {
				std::string val = indexing_parts[i].substr(1, indexing_parts[i].size() - 2);
				std::string dimname;
				if (!datafile.GetDimensionName(Name, data_name_without_index, i, dimname)) {
					return false;
				}

				size_t ix;
				if (!datafile.GetDimensionIx(Name, dimname, val, &ix)) {
					return false;
				}
				data_indices.push_back(ix);
			} else {
				size_t ix;
				try {
					ix = boost::lexical_cast<size_t>(indexing_parts[i]);
				} catch( boost::bad_lexical_cast const& ) {
					LOGERROR("Error converting index to a number");
					return false;
				}
				data_indices.push_back(ix);
			}
		}
	}

	return true;
}
