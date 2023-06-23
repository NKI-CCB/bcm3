#include "Utils.h"
#include "LikelihoodIncucytePopulation.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#define DIRICHLET_APOPTOSIS 1
#define SCURVE_DRUG_EFFECT 1

LikelihoodIncucytePopulation::Well::Well(size_t num_timepoints, size_t num_replicates)
{
	cell_count = VectorReal::Zero(num_timepoints);
	apoptotic_cell_count = VectorReal::Zero(num_timepoints);
	debris = VectorReal::Zero(num_timepoints);
	confluence = MatrixReal::Zero(num_timepoints, num_replicates);
	apoptosis_marker = MatrixReal::Zero(num_timepoints, num_replicates);
}

LikelihoodIncucytePopulation::LikelihoodIncucytePopulation(size_t numthreads, size_t evaluation_threads)
	: evaluation_threads(1)
	, use_pao_control(true)
	, num_combo_drugs(0)
{
	solvers.resize(numthreads);
	parallel_data.resize(numthreads);
}

LikelihoodIncucytePopulation::~LikelihoodIncucytePopulation()
{
}

bool LikelihoodIncucytePopulation::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	drug_name = likelihood_node.get<std::string>("<xmlattr>.drug");
	cell_line = likelihood_node.get<std::string>("<xmlattr>.cell_line");

	bool result = true;

	bcm3::NetCDFDataFile data;
	if (!data.Open("drug_response_data.nc", false)) {
		return false;
	}

	// Find out how many experiments there are for this drug+cell line combination
	std::vector<std::string> data_groups;
	std::string group_name_base = drug_name + std::string("/") + cell_line;
	result &= data.GetSubgroups(group_name_base, data_groups);
	size_t num_experiments = data_groups.size();

	if (num_experiments == 0) {
		LOGERROR("Could not find any experiments for drug \"%s\" and cell line \"%s\" in the data file", drug_name.c_str(), cell_line.c_str());
		return false;
	}

	experiments.resize(num_experiments);
	size_t max_timepoints = 0;
	for (size_t ei = 0; ei < experiments.size(); ei++) {
		Experiment& e = experiments[ei];
		e.experiment_ix = ei;
		e.modeled_cell_treatments.resize(evaluation_threads);

		std::string group = drug_name + std::string("/") + cell_line + std::string("/experiment") + std::to_string((uint64)ei + 1);

		result &= data.GetDimensionSize(group, "time", &e.num_timepoints);
		result &= data.GetDimensionSize(group, "replicates", &e.num_replicates);
		result &= data.GetDimensionSize(group, "drug_concentrations", &e.num_concentrations);

		e.observed_timepoints.resize(e.num_timepoints);
		data.GetValues(group, "time", 0, e.num_timepoints, e.observed_timepoints);

		e.concentrations.resize(e.num_concentrations);
		data.GetValues(group, "drug_concentrations", 0, e.num_concentrations, e.concentrations);
		for (size_t ci = 0; ci < e.num_concentrations; ci++) {
			e.concentrations(ci) = log10(e.concentrations(ci));
		}

		// Load observed data
		e.observed_cell_treatment.drug_treatment.resize(e.num_concentrations, Well(e.num_timepoints, e.num_replicates));
		for (size_t ci = 0; ci < e.num_concentrations; ci++) {
			Well& well = e.observed_cell_treatment.drug_treatment[ci];
			for (size_t ti = 0; ti < e.num_timepoints; ti++) {
				for (size_t repi = 0; repi < e.num_replicates; repi++) {
					data.GetValue(group, "drug_confluence", ti, ci, repi, &well.confluence(ti, repi));
					data.GetValue(group, "drug_apoptosis_marker", ti, ci, repi, &well.apoptosis_marker(ti, repi));
				}
			}
		}
		e.observed_cell_treatment.positive_control = Well(e.num_timepoints, e.num_replicates);
		for (size_t ti = 0; ti < e.num_timepoints; ti++) {
			for (size_t repi = 0; repi < e.num_replicates; repi++) {
				data.GetValue(group, "positive_control_confluence", ti, repi, &e.observed_cell_treatment.positive_control.confluence(ti, repi));
				data.GetValue(group, "positive_control_apoptosis_marker", ti, repi, &e.observed_cell_treatment.positive_control.apoptosis_marker(ti, repi));
			}
		}
		e.observed_cell_treatment.negative_control = Well(e.num_timepoints, e.num_replicates);
		for (size_t ti = 0; ti < e.num_timepoints; ti++) {
			for (size_t repi = 0; repi < e.num_replicates; repi++) {
				data.GetValue(group, "negative_control_confluence", ti, repi, &e.observed_cell_treatment.negative_control.confluence(ti, repi));
				data.GetValue(group, "negative_control_apoptosis_marker", ti, repi, &e.observed_cell_treatment.negative_control.apoptosis_marker(ti, repi));
			}
		}
		e.observed_cell_treatment.ctb.resize(e.num_concentrations);
		for (size_t ci = 0; ci < e.num_concentrations; ci++) {
			data.GetValue(group, "cell_titer_blue_norm", ci, &e.observed_cell_treatment.ctb(ci));
		}

		result &= data.GetAttribute(group, "treatment_time", e.time_of_drug_treatment);
		result &= data.GetAttribute(group, "ctb_time", e.time_of_ctb);
		result &= data.GetAttribute(group, "seeding_density", e.seeding_density);

		// Allocate modeled data structures
		for (size_t threadix = 0; threadix < evaluation_threads; threadix++) {
			e.modeled_cell_treatments[threadix].drug_treatment.resize(e.num_concentrations, Well(e.num_timepoints, 1));
			e.modeled_cell_treatments[threadix].positive_control = Well(e.num_timepoints, 1);
			e.modeled_cell_treatments[threadix].negative_control = Well(e.num_timepoints, 1);
			e.modeled_cell_treatments[threadix].ctb.resize(e.num_concentrations);
		}
		max_timepoints = std::max(max_timepoints, e.num_timepoints);
	}

	// Allocate structures for parallel evaluation
	for (size_t threadix = 0; threadix < evaluation_threads; threadix++) {
		CVODESolverDelay::TDeriviativeFunction derivative = boost::bind(&LikelihoodIncucytePopulation::CalculateDerivative, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4, boost::placeholders::_5, boost::placeholders::_6, boost::placeholders::_7);
		CVODESolverDelay::TJacobianFunction jacobian = boost::bind(&LikelihoodIncucytePopulation::CalculateJacobian, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4, boost::placeholders::_5, boost::placeholders::_6, boost::placeholders::_7, boost::placeholders::_8);
		solvers[threadix].SetDerivativeFunction(derivative);
		solvers[threadix].SetJacobianFunction(jacobian);
		solvers[threadix].Initialize(3, NULL);
		solvers[threadix].SetTolerance(1e-6, 1e-2);

		parallel_data[threadix].initial_conditions.setZero(3);
		parallel_data[threadix].discontinuities_time.setZero(2);
		parallel_data[threadix].output.setZero(3, max_timepoints);
		parallel_data[threadix].delay_y = VectorReal::Zero(3);

		// TODO - random number seed
		parallel_data[threadix].rng.Seed(threadix);
	}

	data.Close();

	return result;
}

bool LikelihoodIncucytePopulation::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	threadix = 0;

	const Real sigma_confluence = values[varset->GetVariableIndex("sigma_confluence")];
	const Real sigma_apoptosis_marker = values[varset->GetVariableIndex("sigma_apoptosis_marker")];
	const Real sigma_ctb = values[varset->GetVariableIndex("sigma_ctb")];

	bool result = true;
	logp = 0.0;

	for (size_t ei = 0; ei < experiments.size(); ei++) {
		Experiment& e = experiments[ei];
		CellTreatment& ct = e.modeled_cell_treatments[threadix];

		result &= SimulateWell(threadix, e, ct.negative_control, false, false, false, std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN(), values);

		if (use_pao_control) {
			const Real pao_proliferation_rate = 0.0;//values[varset->GetVariableIndex("pao_proliferation_rate")];
			const Real pao_apoptosis_rate = values[varset->GetVariableIndex("pao_apoptosis_rate")];
			result &= SimulateWell(threadix, e, ct.positive_control, true, false, false, pao_proliferation_rate, pao_apoptosis_rate, values);
		}

#if SCURVE_DRUG_EFFECT
		for (size_t ci = 0; ci < 9; ci++) {
			Real ic50 = values[varset->GetVariableIndex("drug_proliferation_rate_ic50")];
			Real steepness = exp(values[varset->GetVariableIndex("drug_proliferation_rate_steepness")]);
			Real maxeffect = values[varset->GetVariableIndex("drug_proliferation_rate_maxeffect")];
			Real relative_drug_proliferation_rate = (1.0 - maxeffect) + maxeffect / (exp(steepness * (e.concentrations[ci] - ic50)) + 1.0);
			Real drug_proliferation_rate = values[varset->GetVariableIndex("proliferation_rate")] * relative_drug_proliferation_rate;

			ic50 = values[varset->GetVariableIndex("drug_apoptosis_rate_ic50")];
			steepness = exp(values[varset->GetVariableIndex("drug_apoptosis_rate_steepness")]);
			maxeffect = values[varset->GetVariableIndex("drug_apoptosis_rate_maxeffect")];
			Real drug_apoptosis_rate = maxeffect - maxeffect / (exp(steepness * (e.concentrations[ci] - ic50)) + 1.0);

			result &= SimulateWell(threadix, e, ct.drug_treatment[ci], false, true, false, drug_proliferation_rate, drug_apoptosis_rate, values);
		}
#else
#if DIRICHLET_APOPTOSIS
		Real relative_drug_proliferation_rate = 1.0;
		Real base_apoptosis_rate = values[varset->GetVariableIndex("apoptosis_rate")] * values[varset->GetVariableIndex("proliferation_rate")];
		Real relative_drug_apoptosis_rate = 0.0;
#else
		Real relative_drug_proliferation_rate = 1.0;
		Real drug_apoptosis_rate = values[varset->GetVariableIndex("apoptosis_rate")] * values[varset->GetVariableIndex("proliferation_rate")];
#endif

		//size_t ci = e.num_concentrations - 1;
		size_t ci;
		if (num_combo_drugs == 0) {
			ci = e.num_concentrations - 1;
		} else if (num_combo_drugs == 1) {
			ci = e.num_concentrations - 2;
		} else if (num_combo_drugs == 2) {
			ci = e.num_concentrations - 4;
		} else {
			ASSERT(false);
			LOGERROR("Unsupported number of combination drugs");
			return false;
		}
		while (1) {
			relative_drug_proliferation_rate = (std::max)(relative_drug_proliferation_rate - values[varset->GetVariableIndex(std::string("drug_proliferation_rate_") + std::to_string((uint64)ci + 1))], 0.0);
			Real drug_proliferation_rate = relative_drug_proliferation_rate * values[varset->GetVariableIndex("proliferation_rate")];
#if DIRICHLET_APOPTOSIS
			relative_drug_apoptosis_rate = relative_drug_apoptosis_rate + values[varset->GetVariableIndex(std::string("drug_apoptosis_rate_") + std::to_string((uint64)ci + 1))];
			Real drug_apoptosis_rate = values[varset->GetVariableIndex("drug_apoptosis_rate")] * relative_drug_apoptosis_rate + base_apoptosis_rate;
#else
			drug_apoptosis_rate += values[varset->GetVariableIndex(std::string("drug_apoptosis_rate_") + std::to_string((uint64)ci + 1))];
#endif

			result &= SimulateWell(threadix, e, ct.drug_treatment[ci], false, true, false, drug_proliferation_rate, drug_apoptosis_rate, values);

			if (ci-- == 0) {
				break;
			}
		}
#endif

		if (num_combo_drugs == 1) {
			Real single_drug_proliferation_rate = values[varset->GetVariableIndex("drug_proliferation_rate_single")] * values[varset->GetVariableIndex("proliferation_rate")];
			Real single_drug_apoptosis_rate = values[varset->GetVariableIndex("drug_apoptosis_rate_single")] +
				values[varset->GetVariableIndex("apoptosis_rate")] * values[varset->GetVariableIndex("proliferation_rate")];
			result &= SimulateWell(threadix, e, ct.drug_treatment[e.num_concentrations - 1], false, true, true, single_drug_proliferation_rate, single_drug_apoptosis_rate, values);
		} else if (num_combo_drugs == 2) {
			Real single_drug_proliferation_rate = values[varset->GetVariableIndex("drug_proliferation_rate_single1")] * values[varset->GetVariableIndex("proliferation_rate")];
			Real single_drug_apoptosis_rate = values[varset->GetVariableIndex("drug_apoptosis_rate_single1")] +
				values[varset->GetVariableIndex("apoptosis_rate")] * values[varset->GetVariableIndex("proliferation_rate")];
			result &= SimulateWell(threadix, e, ct.drug_treatment[e.num_concentrations - 1], false, true, true, single_drug_proliferation_rate, single_drug_apoptosis_rate, values);

			single_drug_proliferation_rate = values[varset->GetVariableIndex("drug_proliferation_rate_single2")] * values[varset->GetVariableIndex("proliferation_rate")];
			single_drug_apoptosis_rate = values[varset->GetVariableIndex("drug_apoptosis_rate_single2")] +
				values[varset->GetVariableIndex("apoptosis_rate")] * values[varset->GetVariableIndex("proliferation_rate")];
			result &= SimulateWell(threadix, e, ct.drug_treatment[e.num_concentrations - 2], false, true, true, single_drug_proliferation_rate, single_drug_apoptosis_rate, values);

			single_drug_proliferation_rate = values[varset->GetVariableIndex("drug_proliferation_rate_dual")] * values[varset->GetVariableIndex("proliferation_rate")];
			single_drug_apoptosis_rate = values[varset->GetVariableIndex("drug_apoptosis_rate_dual")] +
				values[varset->GetVariableIndex("apoptosis_rate")] * values[varset->GetVariableIndex("proliferation_rate")];
			result &= SimulateWell(threadix, e, ct.drug_treatment[e.num_concentrations - 3], false, true, true, single_drug_proliferation_rate, single_drug_apoptosis_rate, values);
		}

		// Calculate CTB
		for (size_t ci = 0; ci < e.num_concentrations; ci++) {
			// Assume all cells have died in the positive control
			if (ct.negative_control.cell_count[e.num_timepoints - 1] > 0.0) {
				ct.ctb(ci) = ct.drug_treatment[ci].cell_count[e.num_timepoints - 1] /
					ct.negative_control.cell_count[e.num_timepoints - 1];
			} else {
				// All cells have died in the negative control as well with these parameters - can't calculate CTB
				ct.ctb(ci) = 0.0;
			}
		}

		if (result) {
			// Count the entire time series of all four technical replicates as 1 data point.
			const Real factor = 0.25 / (Real)e.observed_timepoints.size();

			if (use_pao_control) {
				Well& observed_pcw = e.observed_cell_treatment.positive_control;
				Well& modeled_pcw = ct.positive_control;
				for (size_t ti = 0; ti < e.num_timepoints; ti++) {
					for (size_t repi = 0; repi < 4; repi++) {
						logp += factor * bcm3::LogPdfTnu3(observed_pcw.confluence(ti, repi), modeled_pcw.confluence(ti), sigma_confluence, true);
						logp += factor * bcm3::LogPdfTnu3(observed_pcw.apoptosis_marker(ti, repi), modeled_pcw.apoptosis_marker(ti), sigma_apoptosis_marker, true);
					}
				}
			}

			Well& observed_ncw = e.observed_cell_treatment.negative_control;
			Well& modeled_ncw = ct.negative_control;
			for (size_t ti = 0; ti < e.num_timepoints; ti++) {
				for (size_t repi = 0; repi < 4; repi++) {
					logp += factor * bcm3::LogPdfTnu3(observed_ncw.confluence(ti, repi), modeled_ncw.confluence(ti), sigma_confluence, true);
					logp += factor * bcm3::LogPdfTnu3(observed_ncw.apoptosis_marker(ti, repi), modeled_ncw.apoptosis_marker(ti), sigma_apoptosis_marker, true);
				}
			}

			for (size_t ci = 0; ci < e.num_concentrations; ci++) {
				Well& observed_w = e.observed_cell_treatment.drug_treatment[ci];
				Well& modeled_w = ct.drug_treatment[ci];
				for (size_t ti = 0; ti < e.num_timepoints; ti++) {
					for (size_t repi = 0; repi < 4; repi++) {
						logp += factor * bcm3::LogPdfTnu3(observed_w.confluence(ti, repi), modeled_w.confluence(ti), sigma_confluence, true);
						logp += factor * bcm3::LogPdfTnu3(observed_w.apoptosis_marker(ti, repi), modeled_w.apoptosis_marker(ti), sigma_apoptosis_marker, true);
					}
				}

				logp += bcm3::LogPdfTnu3(e.observed_cell_treatment.ctb(ci), ct.ctb(ci), sigma_ctb, true);
			}
		} else {
			logp = -std::numeric_limits<Real>::infinity();
			break;
		}
	}

	return true;
}

bool LikelihoodIncucytePopulation::SimulateWell(size_t threadix, Experiment& e, Well& w, bool pao, bool drug, bool single_drug, Real drug_proliferation_rate, Real drug_apoptosis_rate, const VectorReal& values)
{
	ParallelData& pd = parallel_data[threadix];

	// Cell size is in um; a well is 10.9 mm^2 = 1.09e7 um^2, so cell size as fraction of confluence = size in (um^2 / 1.09e+7) * 100
	pd.cell_size = bcm3::fastpow10(values[varset->GetVariableIndex("log10_cell_size")]) * 9.174312e-6;
	if (pao) {
		pd.apoptotic_cell_size = values[varset->GetVariableIndex("pao_apoptotic_cell_size")] * pd.cell_size;
	} else {
		pd.apoptotic_cell_size = values[varset->GetVariableIndex("apoptotic_cell_size")] * pd.cell_size;
	}
	pd.debris_size = values[varset->GetVariableIndex("debris_size")] * pd.cell_size;
	const Real apoptosis_marker_size = values[varset->GetVariableIndex("apoptosis_marker_size")] * pd.cell_size;
	const Real pao_apoptosis_marker_size = values[varset->GetVariableIndex("pao_apoptosis_marker_size")] * pd.cell_size;
	const Real debris_apoptosis_marker_size = values[varset->GetVariableIndex("debris_apoptosis_marker_size")] * apoptosis_marker_size;

	pd.proliferation_rate = values[varset->GetVariableIndex("proliferation_rate")];
	pd.apoptosis_rate = values[varset->GetVariableIndex("apoptosis_rate")] * pd.proliferation_rate;
	pd.apoptosis_duration = values[varset->GetVariableIndex("apoptosis_duration")];
	pd.apoptosis_remove_rate = values[varset->GetVariableIndex("apoptosis_remove_rate")];
	if (pao) {
		pd.drug_start_time = e.time_of_drug_treatment + values[varset->GetVariableIndex("pao_delay")];;
		pd.drug_effect_time = values[varset->GetVariableIndex("pao_effect_time")];
	} else if (single_drug) {
		pd.drug_start_time = e.time_of_drug_treatment + values[varset->GetVariableIndex("single_drug_delay")];;
		pd.drug_effect_time = values[varset->GetVariableIndex("single_drug_effect_time")];
	} else {
		pd.drug_start_time = e.time_of_drug_treatment + values[varset->GetVariableIndex("drug_delay")];;
		pd.drug_effect_time = values[varset->GetVariableIndex("drug_effect_time")];
	}

	pd.contact_inhibition_start = values[varset->GetVariableIndex("contact_inhibition_start")];
	pd.contact_inhibition_max_confluence = values[varset->GetVariableIndex("contact_inhibition_max_confluence")];
	pd.contact_inhibition_apoptosis_rate = values[varset->GetVariableIndex("contact_inhibition_apoptosis_rate")];

	pd.drug_proliferation_rate = drug_proliferation_rate;
	pd.drug_apoptosis_rate = drug_apoptosis_rate;
	pd.pao = pao;
	pd.drug = drug;

	Real seeding_density_deviation = values[varset->GetVariableIndex(std::string("seeding_density_deviation_") + std::to_string((uint64)e.experiment_ix + 1))];
	Real dead_cell_fraction = values[varset->GetVariableIndex("starting_dead_cell_fraction")];
	
	VectorReal& initial_conditions = pd.initial_conditions;
	initial_conditions(0) = e.seeding_density * bcm3::fastpow10(seeding_density_deviation);
	initial_conditions(1) = dead_cell_fraction * initial_conditions(0);
	initial_conditions(0) -= initial_conditions(1);
	
	VectorReal& discontinuities_time = pd.discontinuities_time;
	if (pao || drug) {
		discontinuities_time.resize(2);
		discontinuities_time(0) = pd.drug_start_time;
		discontinuities_time(1) = pd.drug_start_time + pd.drug_effect_time;
	} else {
		discontinuities_time.resize(0);
	}

	MatrixReal& output = pd.output;
	
	solvers[threadix].SetUserData((void*)threadix);
	bool result = solvers[threadix].Simulate(initial_conditions.data(), e.observed_timepoints, discontinuities_time, output);
	//LOG("Success: %d; steps: %zu", (int)result, solvers[threadix].GetNumStepsTaken());
	if (result) {
		const Real cell_preadherence_size = values[varset->GetVariableIndex("cell_preadherence_size")];
		const Real cell_size_decrease_time = values[varset->GetVariableIndex("cell_adherence_time")];

#if 0
		Real cell_treatment_size_factor, cell_treatment_size_time;
		if (pao) {
			cell_treatment_size_factor = values[varset->GetVariableIndex("pao_size_factor")];
			cell_treatment_size_time = values[varset->GetVariableIndex("pao_size_time")];
		} else if (drug) {
			cell_treatment_size_factor = values[varset->GetVariableIndex("drug_size_factor")];
			cell_treatment_size_time = values[varset->GetVariableIndex("drug_size_time")];
		}
#endif

		for (size_t ti = 0; ti < e.num_timepoints; ti++) {
			w.cell_count[ti] = output(0, ti);
			w.apoptotic_cell_count[ti] = output(1, ti);
			w.debris[ti] = output(2, ti);

			Real cell_size_factor;
			if (e.observed_timepoints(ti) < cell_size_decrease_time) {
				cell_size_factor = cell_preadherence_size + (1.0 - cell_preadherence_size) * e.observed_timepoints(ti) / cell_size_decrease_time;
			} else {
				cell_size_factor = 1.0;
			}
#if 0
			if (pao || drug) {
				if (e.observed_timepoints(ti) > e.time_of_drug_treatment) {
					if (e.observed_timepoints(ti) > e.time_of_drug_treatment + cell_treatment_size_time) {
						cell_size_factor *= cell_treatment_size_factor;
					} else {
						Real t = (e.observed_timepoints(ti) - e.time_of_drug_treatment) / cell_treatment_size_time;
						cell_size_factor *= (1.0 - t) + t * cell_treatment_size_factor;
					}
				} else {
				}
			}
#endif

			w.confluence(ti) = w.cell_count[ti] * pd.cell_size * cell_size_factor +
							   w.apoptotic_cell_count[ti] * pd.apoptotic_cell_size + 
							   w.debris[ti] * pd.debris_size;

			// Apoptosis reagent has been added at same time as drug
			if (e.observed_timepoints(ti) < e.time_of_drug_treatment) {
				w.apoptosis_marker(ti) = 0;
			} else {
				if (pao) {
					w.apoptosis_marker(ti) = w.apoptotic_cell_count[ti] * pao_apoptosis_marker_size +
						w.debris[ti] * debris_apoptosis_marker_size;
				} else {
					w.apoptosis_marker(ti) = w.apoptotic_cell_count[ti] * apoptosis_marker_size +
											 w.debris[ti] * debris_apoptosis_marker_size;
				}
			}
		}
	}
	return result;
}

bool LikelihoodIncucytePopulation::CalculateDerivative(OdeReal t, const OdeReal* y, const std::vector< OdeReal >& history_t, const OdeMatrixReal& history_y, size_t current_dci, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	OdeReal proliferation_rate = (OdeReal)pd.proliferation_rate;
	OdeReal apoptosis_rate = (OdeReal)pd.apoptosis_rate;
	OdeReal apoptosis_remove_rate = (OdeReal)pd.apoptosis_remove_rate;
	if (!CVODESolverDelay::InterpolateHistory(t - (OdeReal)pd.apoptosis_duration, history_t, history_y, pd.delay_y)) {
		return false;
	}
	CalculateDrugEffect(proliferation_rate, apoptosis_rate, y, pd, t, current_dci);
	CalculateContactInhibition(proliferation_rate, y, pd);

	dydt[0] = (proliferation_rate - apoptosis_rate) * y[0];
	dydt[1] = apoptosis_rate * y[0] - apoptosis_remove_rate * pd.delay_y(1);
	dydt[2] = apoptosis_remove_rate * pd.delay_y(1);

	return true;
}

bool LikelihoodIncucytePopulation::CalculateJacobian(OdeReal t, const OdeReal* y, const OdeReal* dydt, const std::vector< OdeReal >& history_t, const OdeMatrixReal& history_y, size_t current_dci, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	OdeReal proliferation_rate = (OdeReal)pd.proliferation_rate;
	OdeReal apoptosis_rate = (OdeReal)pd.apoptosis_rate;
	OdeReal apoptosis_remove_rate = (OdeReal)pd.apoptosis_remove_rate;
#if 0
	if (!CVODESolverDelay::InterpolateHistory(t - (OdeReal)pd.apoptosis_duration, history_t, history_y, pd.delay_y)) {
		return false;
	}
#endif
	CalculateDrugEffect(proliferation_rate, apoptosis_rate, y, pd, t, current_dci);
	CalculateContactInhibition(proliferation_rate, y, pd);

	jac(0, 0) = (proliferation_rate - apoptosis_rate);
	jac(0, 1) = 0.0;
	jac(0, 2) = 0.0;

	jac(1, 0) = apoptosis_rate;
	jac(1, 1) = 0.0;// -apoptosis_remove_rate;
	jac(1, 2) = 0.0;

	jac(2, 0) = 0.0;
	jac(2, 1) = 0.0;//apoptosis_remove_rate;
	jac(2, 2) = 0.0;

	return true;
}

void LikelihoodIncucytePopulation::CalculateDrugEffect(Real& proliferation_rate, Real& apoptosis_rate, const Real* y, const ParallelData& pd, Real t, size_t current_dci)
{
	if (pd.drug || pd.pao) {
		if (current_dci == 0) {
			// No drug yet
		} else if (current_dci == 1) {
			// Increasing drug effect
			OdeReal drug_time_effect = (t - (OdeReal)pd.drug_start_time) / (OdeReal)pd.drug_effect_time;
			proliferation_rate = ((OdeReal)1.0 - drug_time_effect) * proliferation_rate + drug_time_effect * (OdeReal)pd.drug_proliferation_rate;
			apoptosis_rate = ((OdeReal)1.0 - drug_time_effect) * apoptosis_rate + drug_time_effect * (OdeReal)pd.drug_apoptosis_rate;
		} else if (current_dci == 2) {
			// Full effect
			proliferation_rate = (OdeReal)pd.drug_proliferation_rate;
			apoptosis_rate = (OdeReal)pd.drug_apoptosis_rate;
		}
	}
}

void LikelihoodIncucytePopulation::CalculateContactInhibition(Real& proliferation_rate, const Real* y, const ParallelData& pd)
{
	OdeReal confluence = (OdeReal)(0.01 * (y[0] * pd.cell_size + y[1] * pd.apoptotic_cell_size + y[2] * pd.debris_size));
	if (confluence > pd.contact_inhibition_start) {
		// ci = strength of contact inhbition, ci == 1 means the growth is fully inhibited.
		OdeReal ci = (OdeReal)((confluence - pd.contact_inhibition_start) / (pd.contact_inhibition_max_confluence - pd.contact_inhibition_start));
		ci = (std::max)((std::min)(ci, 1.0), 0.0);
		//OdeReal proliferation_ci = ci * (OdeReal)effective_proliferation_rate;
		//effective_proliferation_rate -= proliferation_ci * (OdeReal)(1.0 - pd.contact_inhibition_apoptosis_rate);
		//effective_apoptosis_rate += proliferation_ci * (OdeReal)pd.contact_inhibition_apoptosis_rate;
		proliferation_rate *= (1.0 - ci);
	}
}

MatrixReal LikelihoodIncucytePopulation::GetSimulatedCellCount(size_t experiment_ix) const
{
	const Experiment& e = experiments[experiment_ix];
	MatrixReal trajectories(e.num_timepoints, e.num_concentrations + 2);

	for (ptrdiff_t i = 0; i < e.num_timepoints; i++) {
		for (ptrdiff_t j = 0; j < e.num_concentrations; j++) {
			trajectories.col(j) = e.modeled_cell_treatments[0].drug_treatment[j].cell_count;
		}
		trajectories.col(e.num_concentrations + 0) = e.modeled_cell_treatments[0].positive_control.cell_count;
		trajectories.col(e.num_concentrations + 1) = e.modeled_cell_treatments[0].negative_control.cell_count;
	}

	return trajectories;
}

MatrixReal LikelihoodIncucytePopulation::GetSimulatedApoptoticCellCount(size_t experiment_ix) const
{
	const Experiment& e = experiments[experiment_ix];
	MatrixReal trajectories(e.num_timepoints, e.num_concentrations + 2);

	for (ptrdiff_t i = 0; i < e.num_timepoints; i++) {
		for (ptrdiff_t j = 0; j < e.num_concentrations; j++) {
			trajectories.col(j) = e.modeled_cell_treatments[0].drug_treatment[j].apoptotic_cell_count;
		}
		trajectories.col(e.num_concentrations + 0) = e.modeled_cell_treatments[0].positive_control.apoptotic_cell_count;
		trajectories.col(e.num_concentrations + 1) = e.modeled_cell_treatments[0].negative_control.apoptotic_cell_count;
	}

	return trajectories;
}

MatrixReal LikelihoodIncucytePopulation::GetSimulatedDebris(size_t experiment_ix) const
{
	const Experiment& e = experiments[experiment_ix];
	MatrixReal trajectories(e.num_timepoints, e.num_concentrations + 2);

	for (ptrdiff_t i = 0; i < e.num_timepoints; i++) {
		for (ptrdiff_t j = 0; j < e.num_concentrations; j++) {
			trajectories.col(j) = e.modeled_cell_treatments[0].drug_treatment[j].debris;
		}
		trajectories.col(e.num_concentrations + 0) = e.modeled_cell_treatments[0].positive_control.debris;
		trajectories.col(e.num_concentrations + 1) = e.modeled_cell_treatments[0].negative_control.debris;
	}

	return trajectories;
}

MatrixReal LikelihoodIncucytePopulation::GetSimulatedConfluence(size_t experiment_ix) const
{
	const Experiment& e = experiments[experiment_ix];
	MatrixReal trajectories(e.num_timepoints, e.num_concentrations + 2);

	for (ptrdiff_t i = 0; i < e.num_timepoints; i++) {
		for (ptrdiff_t j = 0; j < e.num_concentrations; j++) {
			trajectories.col(j) = e.modeled_cell_treatments[0].drug_treatment[j].confluence;
		}
		trajectories.col(e.num_concentrations + 0) = e.modeled_cell_treatments[0].positive_control.confluence;
		trajectories.col(e.num_concentrations + 1) = e.modeled_cell_treatments[0].negative_control.confluence;
	}

	return trajectories;
}

MatrixReal LikelihoodIncucytePopulation::GetSimulatedApoptosisMarker(size_t experiment_ix) const
{
	const Experiment& e = experiments[experiment_ix];
	MatrixReal trajectories(e.num_timepoints, e.num_concentrations + 2);

	for (ptrdiff_t i = 0; i < e.num_timepoints; i++) {
		for (ptrdiff_t j = 0; j < e.num_concentrations; j++) {
			trajectories.col(j) = e.modeled_cell_treatments[0].drug_treatment[j].apoptosis_marker;
		}
		trajectories.col(e.num_concentrations + 0) = e.modeled_cell_treatments[0].positive_control.apoptosis_marker;
		trajectories.col(e.num_concentrations + 1) = e.modeled_cell_treatments[0].negative_control.apoptosis_marker;
	}

	return trajectories;
}

VectorReal LikelihoodIncucytePopulation::GetSimulatedCTB(size_t experiment_ix) const
{
	const Experiment& e = experiments[experiment_ix];
	return e.modeled_cell_treatments[0].ctb;
}
