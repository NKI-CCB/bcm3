#include "Utils.h"
#include "fISAExperimentDrugRange.h"
#include "fISAExperimentSingleCondition.h"
#include "SignalingNetwork.h"
#include "ProbabilityDistributions.h"

#if TODO

fISAExperimentDrugRange::fISAExperimentDrugRange(size_t numthreads)
	: fISAExperiment(numthreads)
{
}

fISAExperimentDrugRange::~fISAExperimentDrugRange()
{
}

bool fISAExperimentDrugRange::LoadTypeSpecificNodes(const boost::property_tree::ptree& xml_node, const bcm3::NetCDFDataFile& datafile)
{
	bool result = true;

	try {
		const boost::property_tree::ptree& drug_range = xml_node.get_child("drug_range");

		if (!drug_species_name.empty()) {
			LOGERROR("There can be only one drug_range data in an experiment");
			return false;
		}

		drug_species_name = drug_range.get<std::string>("<xmlattr>.species_name");
		drug_model_ix = network->GetSignalingMoleculeIxByName(drug_species_name);
		if (drug_model_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Could not find signaling node \"%s\" for drug range data", drug_species_name.c_str());
			return false;
		}

		std::string data_name = drug_range.get<std::string>("<xmlattr>.concentrations_data_name");
				
		size_t concentration_count = 0;
		result &= datafile.GetDimensionSize(Name, data_name, &concentration_count);
		result &= datafile.GetValues(Name, data_name, 0, concentration_count, drug_concentrations);
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return result;
}

bool fISAExperimentDrugRange::InitializeParallelData()
{
	parallel_data.resize(numthreads);
	for (size_t threadix = 0; threadix < numthreads; threadix++) {
		parallel_data[threadix].transformed_values.resize(varset->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());
		parallel_data[threadix].activities.resize(cell_lines.size());
		parallel_data[threadix].modeled_data_values.resize(cell_lines.size());
		for (size_t ci = 0; ci < cell_lines.size(); ci++) {
			parallel_data[threadix].activities[ci].resize(drug_concentrations.size(), VectorReal::Constant(network->GetMoleculeCount(), std::numeric_limits<Real>::quiet_NaN()));
			parallel_data[threadix].modeled_data_values[ci].resize(drug_concentrations.size(), VectorReal::Constant(data.size(), std::numeric_limits<Real>::quiet_NaN()));
		}
	}

	return true;
}

bool fISAExperimentDrugRange::ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile)
{
	bool result = true;

	data.push_back(DataPart());
	DataPart& p = *data.rbegin();
	
	if (!ParseDataPartBase(xml_node, p, other_experiments)) {
		return false;
	}
	
	std::string data_name_without_index;
	std::vector<size_t> data_indices;
	if (!ParseDataFileReference(p.data_name, data_name_without_index, data_indices, datafile)) {
		LOGERROR("Unable to parse data reference \"%s\" for data node", p.data_name.c_str());
		return false;
	}
	
	size_t data_dim_count = 0;
	datafile.GetDimensionCount(Name, data_name_without_index, &data_dim_count);
	if (data_dim_count != data_indices.size()) {
		LOGERROR("The data reference \"%s\" in the likelihood has %u dimensions, but the data entry in the data file has %u dimensions", p.data_name.c_str(), data_indices.size(), data_dim_count);
		return false;
	}
	
	size_t cell_line_count;
	std::string cell_line_dim_name;
	result &= datafile.GetDimensionName(Name, data_name_without_index, 1, cell_line_dim_name);
	result &= datafile.GetDimensionSize(Name, cell_line_dim_name, &cell_line_count);
	if (cell_line_count != cell_lines.size()) {
		LOGERROR("Inconsistent number of cell lines for data variable \"%s\"", p.data_name.c_str());
		return false;
	}

	size_t drug_concentration_count;
	std::string drug_concentration_dim_name;
	result &= datafile.GetDimensionName(Name, data_name_without_index, 2, drug_concentration_dim_name);
	result &= datafile.GetDimensionSize(Name, drug_concentration_dim_name, &drug_concentration_count);
	if ((int)drug_concentration_count != drug_concentrations.size()) {
		LOGERROR("Inconsistent number of drug concentrations for data variable \"%s\"", p.data_name.c_str());
		return false;
	}
	
	if (data_indices.size() == 3) {
		// No replicates
		p.data.resize(1);
		p.data[0] = MatrixReal::Constant(cell_lines.size(), drug_concentrations.size(), std::numeric_limits<Real>::quiet_NaN());
		for (int i = 0; i < drug_concentrations.size(); i++) {
			for (size_t ci = 0; ci < cell_lines.size(); ci++) {
				result &= datafile.GetValue(Name, data_name_without_index, data_indices[0], ci, (size_t)i, 0, &p.data[0](ci, i));
			}
		}
	} else if (data_indices.size() == 4) {
		// There are replicates
		size_t num_replicates;
		std::string replicate_dim_name;
		result &= datafile.GetDimensionName(Name, data_name_without_index, 3, replicate_dim_name);
		result &= datafile.GetDimensionSize(Name, replicate_dim_name, &num_replicates);

		p.data.resize(num_replicates);
		for (size_t ri = 0; ri < num_replicates; ri++) {
			p.data[ri] = MatrixReal::Constant(cell_lines.size(), drug_concentrations.size(), std::numeric_limits<Real>::quiet_NaN());
			for (int i = 0; i < drug_concentrations.size(); i++) {
				for (size_t ci = 0; ci < cell_lines.size(); ci++) {
					result &= datafile.GetValue(Name, data_name_without_index, data_indices[0], ci, (size_t)i, ri, &p.data[ri](ci, i));
				}
			}
		}
	} else {
		LOGERROR("Data in a drug range experiment should have 3 or 4 dimensions, but data reference \"%s\" has %u dimensions", p.data_name.c_str(), data_indices.size());
		return false;
	}

	return result;
}

bool fISAExperimentDrugRange::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	ParallelData& pd = parallel_data[threadix];

	for (size_t vi = 0; vi < varset->GetNumVariables(); vi++) {
		pd.transformed_values[vi] = varset->TransformVariable(vi, values[vi]);
	}

	const Real tdist_nu3_C = bcm3::LogPdfT_CalcC(3.0);

	logp = 0.0;

	for (int dci = 0; dci < drug_concentrations.size(); dci++) {
		for (size_t ci = 0; ci < cell_lines.size(); ci++) {
			VectorReal& activities = pd.activities[ci][dci];

			// Calculate model activities
			PrepareActivitiesCalculation(activities, pd.expression, pd.transformed_values.data(), ci);
			activities(drug_model_ix) = drug_concentrations(dci);
			if (network->Calculate(threadix, activities, pd.expression, pd.transformed_values.data())) {
				// Evaluate data likelihood
				for (size_t di = 0; di < data.size(); di++) {
					const DataPart& d = data[di];

					const Real var = d.fixed_sd ? d.fixed_sd_value : pd.transformed_values[d.parameter_sd_ix];

					Real y = activities(d.model_ix);
					if (d.model_ix2 != std::numeric_limits<size_t>::max()) {
						Real sum_parameter = pd.transformed_values[d.data_sum_parameter_ix];
						y = sum_parameter * y + (1.0 - sum_parameter) * activities(d.model_ix2);
					}
					if (d.data_is_inactive_form) {
						y = pd.expression(d.model_ix) - y;
					}
					if (d.use_scale) {
						y *= pd.transformed_values[d.parameter_scale_ix];
					}
					if (d.use_base) {
						if (d.fixed_base) {
							y += d.fixed_base_value;
						}
						else {
							y += pd.transformed_values[d.parameter_base_ix];
						}
					}

					if (!d.relative_to_experiment.empty()) {
						// This data is relative to some other experiment; find the other experiment
						for (std::vector< std::unique_ptr<fISAExperiment> >::const_iterator ei = other_experiments.begin(); ei != other_experiments.end(); ++ei) {
							if ((*ei)->GetName() == d.relative_to_experiment) {
								fISAExperimentSingleCondition* expsc = dynamic_cast<fISAExperimentSingleCondition*>(ei->get());
								ASSERT(expsc != NULL);
								const Real relative_reference_y = expsc->parallel_data[threadix].activities[ci](d.model_ix);

								// For now, this is used only for the proliferation node, in which case we want the function 
								// cell_titer = exp((growth_rate - untreated_growth_rate) * 72 * 0.05)
								// The 72 comes from 3 day drug treatment, 3*24 hours
								// The 0.05 comes from: proliferation is between 0 and 1, we'll set 1 to a growth rate of 0.05 (doubling time of ~14 hours)
								// Should probably make this more generic at some point.
								y = exp((y - relative_reference_y) * 72.0 * 0.05);
								break;
							}
						}
					}

					pd.modeled_data_values[ci][dci](di) = y;
					for (size_t ri = 0; ri < d.data.size(); ri++) {
						logp += d.weight * bcm3::LogPdfT(d.data[ri](ci, dci), y, var, 3.0, tdist_nu3_C, true);
					}
				}
			} else {
				logp = -std::numeric_limits<Real>::infinity();
				break;
			}
		}

		if (logp == -std::numeric_limits<Real>::infinity()) {
			break;
		}
	}

	return true;
}

#endif
