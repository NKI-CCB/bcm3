#include "Utils.h"
#include "fISAExperimentSingleCondition.h"
#include "SignalingNetwork.h"
#include "ProbabilityDistributions.h"
#include "TaskManager.h"

fISAExperimentSingleCondition::fISAExperimentSingleCondition()
{
}

fISAExperimentSingleCondition::~fISAExperimentSingleCondition()
{
}

bool fISAExperimentSingleCondition::StartEvaluateLogProbability(const VectorReal& values, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	for (size_t vi = 0; vi < varset->GetNumVariables(); vi++) {
		transformed_values[vi] = varset->TransformVariable(vi, values[vi]);
		if (transformed_values[vi] == std::numeric_limits<Real>::infinity()) {
			LOGWARNING("Overflow of parameter %u after transformation of value %g; truncating to highest floating point value", vi, values[vi]);
			transformed_values[vi] = std::numeric_limits<Real>::max();
		} else if (transformed_values[vi] == -std::numeric_limits<Real>::infinity()) {
			LOGWARNING("Negative overflow of parameter %u after transformation of value %g; truncating to lowest floating point value", vi, values[vi]);
			transformed_values[vi] = std::numeric_limits<Real>::lowest();
		}
	}

	if (network->GetNumEvaluationThreads() > 1) {
		for (size_t ci = 0; ci < cell_lines.size(); ci++) {
			bcm3::TTask task = boost::bind(&fISAExperimentSingleCondition::EvaluateCellLine, this, boost::placeholders::_1, boost::placeholders::_2);
			evaluation_tasks[ci] = task_manager->AddTask(task, (void*)ci);
		}
	}

	return true;
}

bool fISAExperimentSingleCondition::FinishEvaluateLogProbability(const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	logp = 0.0;
	for (size_t ci = 0; ci < cell_lines.size(); ci++) {
		bool result = true;
		if (network->GetNumEvaluationThreads() > 1) {
			result = task_manager->WaitTask(evaluation_tasks[ci]);
		} else {
			result = fISAExperimentSingleCondition::EvaluateCellLine((void*)ci, 0);
		}
		if (result) {
			logp += cell_line_logp(ci);
		} else {
			LOGERROR("Evaluation %u failed", ci);
			return false;
		}
	}
	
	return true;
}

void fISAExperimentSingleCondition::GetObservedData(size_t data_ix, MatrixReal& out_values) const
{
	const DataPart& d = data[data_ix];

	out_values.resize(cell_lines.size(), d.data.size());
	for (size_t ci = 0; ci < cell_lines.size(); ci++) {
		for (size_t ri = 0; ri < d.data.size(); ri++) {
			out_values(ci,ri) = d.data[ri](ci);
		}
	}
}

void fISAExperimentSingleCondition::GetModeledData(size_t threadix, size_t data_ix, VectorReal& out_values) const
{
	const ParallelData& pd = parallel_data[threadix];

	out_values.resize(cell_lines.size());
	for (size_t ci = 0; ci < cell_lines.size(); ci++) {
		out_values(ci) = stored_modeled_data_values[ci](data_ix);
	}
}

void fISAExperimentSingleCondition::GetModeledActivities(size_t threadix, MatrixReal& out_values) const
{
	const ParallelData& pd = parallel_data[threadix];

	out_values.resize(cell_lines.size(), network->GetMoleculeCount());
	for (size_t ci = 0; ci < cell_lines.size(); ci++) {
		out_values.row(ci) = stored_activities[ci];
	}
}

bool fISAExperimentSingleCondition::InitializeParallelData(size_t evaluation_threads)
{
	parallel_data.resize(evaluation_threads);
	for (size_t threadix = 0; threadix < evaluation_threads; threadix++) {
		parallel_data[threadix].multi_activities.resize(network->GetMultirootSolveCount(), VectorReal::Constant(network->GetMoleculeCount(), std::numeric_limits<Real>::quiet_NaN()));
		parallel_data[threadix].multi_modeled_data_values.resize(network->GetMultirootSolveCount(), VectorReal::Constant(data.size(), std::numeric_limits<Real>::quiet_NaN()));
		parallel_data[threadix].multisolve_logp.setZero(network->GetMultirootSolveCount());
	}

	cell_line_logp.setZero(cell_lines.size());
	evaluation_tasks.resize(cell_lines.size(), 0);
	stored_activities.resize(cell_lines.size());
	stored_modeled_data_values.resize(cell_lines.size());
	for (size_t ci = 0; ci < cell_lines.size(); ci++) {
		stored_activities[ci] = VectorReal::Constant(network->GetMoleculeCount(), std::numeric_limits<Real>::quiet_NaN());
		stored_modeled_data_values[ci] = VectorReal::Constant(data.size(), std::numeric_limits<Real>::quiet_NaN());
	}

	return true;
}

bool fISAExperimentSingleCondition::ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile)
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

	// Assumptions:
	// If experiment is regular:
	// - Dimension 1 is the gene/epitope/etc
	// - Dimension 2 are the cell lines
	// - Dimension 3 are the replicates, if any
	if (data_indices.size() == 2) {
		// No replicates
		p.data.resize(1);
		result &= datafile.GetValuesDim2(Name, data_name_without_index, data_indices[0], 0, cell_lines.size(), p.data[0]);
	} else if (data_indices.size() == 3) {
		// There are replicates
		size_t num_replicates;
		std::string replicate_dim_name;
		result &= datafile.GetDimensionName(Name, data_name_without_index, 2, replicate_dim_name);
		result &= datafile.GetDimensionSize(Name, replicate_dim_name, &num_replicates);

		p.data.resize(num_replicates);
		for (size_t ri = 0; ri < num_replicates; ri++) {
			result &= datafile.GetValuesDim2(Name, data_name_without_index, data_indices[0], 0, ri, cell_lines.size(), p.data[ri]);
		}

		if (!result) {
			return false;
		}
	} else {
		LOGERROR("Data in a single condition experiment should have 2 or 3 dimensions, but data reference \"%s\" has %u dimensions", p.data_name.c_str(), data_indices.size());
		return false;
	}

	return true;
}

bool fISAExperimentSingleCondition::EvaluateCellLine(void* cell_line_ix_as_ptr, size_t eval_thread_ix)
{
	ParallelData& pd = parallel_data[eval_thread_ix];
	size_t cell_line_ix = (size_t)cell_line_ix_as_ptr;

	// Calculate model activities
	PrepareActivitiesCalculation(pd.multi_activities[0], pd.expression, transformed_values.data(), cell_line_ix);
	for (size_t mi = 1; mi < network->GetMultirootSolveCount(); mi++) {
		pd.multi_activities[mi] = pd.multi_activities[0];
	}

	if (!network->Calculate(eval_thread_ix, pd.multi_activities, pd.expression, transformed_values.data())) {
		cell_line_logp(cell_line_ix) = -std::numeric_limits<Real>::infinity();
		return true;
	}
	
	// Evaluate data likelihood
	pd.multisolve_logp.setZero();

	for (size_t mi = 0; mi < network->GetMultirootSolveCount(); mi++) {
//		LOG("Experiment %s - cell line %u - solve %u - MTORC2 activity=%g; AKT_S473 activity=%g, PI3K activity=%g, IRS1 activity=%g", Name.c_str(), cell_line_ix, mi,
//			pd.multi_activities[mi](network->GetSignalingMoleculeIxByName("MTORC2")),
//			pd.multi_activities[mi](network->GetSignalingMoleculeIxByName("AKT_S473")),
//			pd.multi_activities[mi](network->GetSignalingMoleculeIxByName("IRS1")),
//			pd.multi_activities[mi](network->GetSignalingMoleculeIxByName("PI3K")));
		for (size_t di = 0; di < data.size(); di++) {
			const DataPart& d = data[di];
			Real z;
			if (d.model_ix2 != std::numeric_limits<size_t>::max()) {
				Real sum_parameter = 0.5;
				if (d.data_sum_parameter_ix != std::numeric_limits<size_t>::max()) {
					sum_parameter = transformed_values[d.data_sum_parameter_ix];
				}
				z = sum_parameter * pd.multi_activities[mi](d.model_ix) + (1.0 - sum_parameter) * pd.multi_activities[mi](d.model_ix2);
			} else {
				z = pd.multi_activities[mi](d.model_ix);
			}
			if (d.data_is_inactive_form) {
				Real me = network->max_expression_function(pd.expression(d.model_ix), d.model_ix, transformed_values.data());
				z = me - z;
			}
			if (d.expression_ix != std::numeric_limits<size_t>::max()) {
				if (d.expression_ix2 != std::numeric_limits<size_t>::max()) {
					Real sum_parameter = 0.5;
					if (d.expression_sum_parameter_ix != std::numeric_limits<size_t>::max()) {
						sum_parameter = transformed_values[d.expression_sum_parameter_ix];
					}
					z *= sum_parameter * pd.expression(d.expression_ix) + (1.0 - sum_parameter) * pd.expression(d.expression_ix2);
				} else {
					z *= pd.expression(d.expression_ix);
				}
			}
			if (d.use_scale) {
				if (d.scale_per_cell_line) {
					z *= transformed_values[d.parameter_scale_ix + cell_line_ix];
				} else {
					z *= transformed_values[d.parameter_scale_ix];
				}
			}
			if (d.use_base) {
				if (d.fixed_base) {
					z += d.fixed_base_value;
				} else {
					z += transformed_values[d.parameter_base_ix];
				}
			}

			Real thisp = 0.0;
			if (d.likelihood_fn == LF_Normal || d.likelihood_fn == LF_TruncatedNormal || d.likelihood_fn == LF_StudentT || d.likelihood_fn == LF_TruncatedStudentT) {
				Real var = d.fixed_sd ? d.fixed_sd_value : transformed_values[d.parameter_sd_ix];
				if (d.scale_var_with_mean) {
					if (d.parameter_sd_base_ix != std::numeric_limits<size_t>::max()) {
						if (d.fixed_sd_scale_value == d.fixed_sd_scale_value) {
							var = var + d.fixed_sd_scale_value * fabs(z);
						} else {
							Real varbase = transformed_values[d.parameter_sd_base_ix];
							var = varbase + var * fabs(z);
						}
					} else {
						var = var * fabs(z);
					}
				}

				if (d.likelihood_fn == LF_Normal) {
					for (size_t ri = 0; ri < d.data.size(); ri++) {
						thisp += bcm3::LogPdfNormal(d.data[ri](cell_line_ix), z, var, true);
					}
				} else if (d.likelihood_fn == LF_TruncatedNormal) {
					for (size_t ri = 0; ri < d.data.size(); ri++) {
						thisp += bcm3::LogPdfTruncatedNormal(d.data[ri](cell_line_ix), z, var, 0.0, 1.0, true);
					}
				} else if (d.likelihood_fn == LF_TruncatedStudentT) {
					z = std::min(z, 1.0);
					for (size_t ri = 0; ri < d.data.size(); ri++) {
						thisp += bcm3::LogTruncatedPdfTnu3(d.data[ri](cell_line_ix), z, var, 0.0, 1.0, true);
					}
				} else {
					for (size_t ri = 0; ri < d.data.size(); ri++) {
						thisp += bcm3::LogPdfTnu3(d.data[ri](cell_line_ix), z, var, true);
					}
				}
			} else if (d.likelihood_fn == LF_Binomial) {
				for (size_t ri = 0; ri < d.data.size(); ri++) {
					Real x = d.data[ri](cell_line_ix);
					if (x == x) {
						if (x == 0.0) {
							if (z >= 1.0) {
								thisp = -std::numeric_limits<Real>::infinity();
							} else {
								thisp += log(1.0 - z);
							}
						} else if (x == 1.0) {
							if (z <= 0.0) {
								thisp = -std::numeric_limits<Real>::infinity();
							} else {
								thisp += log(z);
							}
						} else {
							LOGERROR("Invalid data value %g for binomial, must be 0 or 1", x);
							return false;
						}
					}
				}
			} else if (d.likelihood_fn == LF_Beta) {
				const Real precision = d.fixed_sd ? d.fixed_precision_value : transformed_values[d.parameter_precision_ix];
				for (size_t ri = 0; ri < d.data.size(); ri++) {
					const Real nu = 3.0;

					Real x = d.data[ri](cell_line_ix);
					if (x == x) {
						if (x == 1.0) {
							x = 0.99;
						} else if (x == 0.0) {
							x = 0.01;
						}
						thisp += bcm3::LogPdfBetaPrecision(x, z, precision);
					}
				}
			} else if (d.likelihood_fn == LF_OrderedProbit || d.likelihood_fn == LF_OrderedRobit) {
				for (size_t ri = 0; ri < d.data.size(); ri++) {
#if 0
					const Real nri_values[] = { 0, 0.12, 0.14, 0.17, 0.2, 0.25, 0.29, 0.33, 0.38, 0.4, 0.43, 0.5, 0.57, 0.6, 0.62, 0.67, 0.71, 0.75, 0.8, 0.83, 0.86, 0.88, 1 };
					Real cumulative_thresholds[22];
					Real prev = 0.0;
					for (size_t i = 0; i < d.ordered_probit_thresholds.size(); i++) {
						cumulative_thresholds[i] = prev + transformed_values[d.ordered_probit_thresholds[i]];
						prev = cumulative_thresholds[i];
						cumulative_thresholds[i] = (i + 1) / (Real)(d.ordered_probit_thresholds.size() + 1);
					}
#else
					//const Real nri_values[] = { 1.0/3.0, 2.0/3.0, 1.0 };
					//const Real nri_values[] = { 0.3, 0.7, 1.0 };
					//const Real nri_values[] = { 0.25, 0.5, 0.75 };
					const Real nri_values[] = { 0.4, 0.6 };
					Real cumulative_thresholds[3];
					Real prev = 0.0;
					//cumulative_thresholds[0] = 0.25;
					for (size_t i = 0; i < d.ordered_probit_thresholds.size(); i++) {
						cumulative_thresholds[i] = prev + transformed_values[d.ordered_probit_thresholds[i]];
						prev = cumulative_thresholds[i];
					}
					//cumulative_thresholds[0] = 0.7095834;
					//cumulative_thresholds[1] = 1.5255129;
					//cumulative_thresholds[2] = 3.6134332;
#endif

					size_t var_ix = varset->GetVariableIndex("nri_ordered_probit_sigma", false);
					Real variance;
					if (var_ix == std::numeric_limits<size_t>::max()) {
						variance = 1.0;
					} else {
						variance = transformed_values[var_ix];
					}
					Real nu;
					if (d.likelihood_fn == LF_OrderedRobit) {
						size_t var_ix_nu = varset->GetVariableIndex("nri_ordered_robit_nu", false);
						if (var_ix_nu == std::numeric_limits<size_t>::max()) {
							nu = 3.0;
						} else {
							nu = transformed_values[var_ix_nu];
						}
					} else {
						nu = 1.0;
					}

					Real x = d.data[ri](cell_line_ix);
					if (x == x) {
						size_t ix = d.ordered_probit_thresholds.size();
						for (size_t i = 0; i < d.ordered_probit_thresholds.size(); i++) {
							if (x < nri_values[i]) {
								//if (x >= nri_values[i] - 1e-5 && x <= nri_values[i] + 1e-5) {
								ix = i;
								break;
							}
						}

						if (d.likelihood_fn == LF_OrderedRobit) {
							if (ix == 0) {
								Real threshold = cumulative_thresholds[ix];
								thisp += log(bcm3::CdfT(threshold, z, variance, nu));
							} else if (ix == d.ordered_probit_thresholds.size()) {
								Real threshold = cumulative_thresholds[ix - 1];
								thisp += log(1 - bcm3::CdfT(threshold, z, variance, nu));
							} else {
								Real threshold1 = cumulative_thresholds[ix];
								Real threshold2 = cumulative_thresholds[ix - 1];
								thisp += log(bcm3::CdfT(threshold1, z, variance, nu) - bcm3::CdfT(threshold2, z, variance, nu));
							}
						} else if (d.likelihood_fn == LF_OrderedProbit) {
							if (ix == 0) {
								Real threshold = cumulative_thresholds[ix];
								Real p = bcm3::CdfNormal(threshold, z, variance);
								thisp += log(p);
							} else if (ix == d.ordered_probit_thresholds.size()) {
								Real threshold = cumulative_thresholds[ix - 1];
								Real p = 1 - bcm3::CdfNormal(threshold, z, variance);
								thisp += log(p);
							} else {
								Real threshold1 = cumulative_thresholds[ix];
								Real threshold2 = cumulative_thresholds[ix - 1];
								Real p = bcm3::CdfNormal(threshold1, z, variance) - bcm3::CdfNormal(threshold2, z, variance);
								thisp += log(p);
							}
						}
					}
				}
			}

			pd.multi_modeled_data_values[mi](di) = z;

			pd.multisolve_logp(mi) += d.weight * thisp;
		}
	}

	int best_ix;
	Real bestp = pd.multisolve_logp.maxCoeff(&best_ix);

	if (bestp == -std::numeric_limits<Real>::infinity()) {
		cell_line_logp(cell_line_ix) = -std::numeric_limits<Real>::infinity();
	} else {
		//LOG("best_ix = %d, %.15g, %.15g", best_ix, cell_line_p(0), bestp);
		//for (int i = 0; i < network->GetMoleculeCount(); i++) {
		//	LOG("%g\t%g", parallel_data[threadix].multi_activities[cell_line_ix][0][i], parallel_data[threadix].multi_activities[cell_line_ix][best_ix][i]);
		//}

		cell_line_logp(cell_line_ix) = bestp;
		stored_activities[cell_line_ix] = pd.multi_activities[best_ix];
		stored_modeled_data_values[cell_line_ix] = pd.multi_modeled_data_values[best_ix];
	}
	return true;
}
