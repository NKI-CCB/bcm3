#include "Utils.h"
#include "LikelihoodPopPKTrajectory.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include "ODESolverDP5.h"

#include <boost/property_tree/xml_parser.hpp>
using namespace boost::placeholders;

typedef std::function<bool(OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;

LikelihoodPopPKTrajectory::ParallelData::ParallelData()
	: dose(std::numeric_limits<Real>::quiet_NaN())
	, dosing_interval(std::numeric_limits<Real>::quiet_NaN())
	, dose_after_dose_change(std::numeric_limits<Real>::quiet_NaN())
	, dose_change_time(std::numeric_limits<Real>::quiet_NaN())
	, intermittent(0)
	, k_absorption(std::numeric_limits<Real>::quiet_NaN())
	, k_excretion(std::numeric_limits<Real>::quiet_NaN())
	, k_elimination(std::numeric_limits<Real>::quiet_NaN())
	, k_vod(std::numeric_limits<Real>::quiet_NaN())
	, k_intercompartmental(std::numeric_limits<Real>::quiet_NaN())
	, k_transit(std::numeric_limits<Real>::quiet_NaN())
	, k_biphasic_switch_time(std::numeric_limits<Real>::quiet_NaN())
	, k_absorption2(std::numeric_limits<Real>::quiet_NaN())
	, biphasic_switch(std::numeric_limits<Real>::quiet_NaN())
	, current_dose_time(std::numeric_limits<Real>::quiet_NaN())
{
}

LikelihoodPopPKTrajectory::LikelihoodPopPKTrajectory(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, pk_type(PKMT_Undefined)
{
	solvers.resize(sampling_threads);
	parallel_data.resize(sampling_threads);
}

LikelihoodPopPKTrajectory::~LikelihoodPopPKTrajectory()
{
}

bool LikelihoodPopPKTrajectory::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	bool result = true;

	this->varset = varset;
	std::string trial;
	std::string pkdata_file;

	try {
		boost::property_tree::ptree& modelnode = likelihood_node.get_child("pk_model");
		drug = modelnode.get<std::string>("<xmlattr>.drug");
		std::string pk_type_str = modelnode.get<std::string>("<xmlattr>.type");
		trial = modelnode.get<std::string>("<xmlattr>.trial");
		pkdata_file = modelnode.get<std::string>("<xmlattr>.pkdata_file");

		if (pk_type_str == "one") {
			pk_type = PKMT_OneCompartment;
		} else if (pk_type_str == "two") {
			pk_type = PKMT_TwoCompartment;
		} else if (pk_type_str == "two_biphasic") {
			pk_type = PKMT_TwoCompartmentBiPhasic;
		} else if (pk_type_str == "one_transit") {
			pk_type = PKMT_OneCompartmentTransit;
		} else if (pk_type_str == "two_transit") {
			pk_type = PKMT_TwoCompartmentTransit;
		} else {
			pk_type = PKMT_Undefined;
		}
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	bcm3::NetCDFDataFile data;
	if (!data.Open(pkdata_file, false)) {
		return false;
	}

	size_t num_timepoints, num_patients;
	result &= data.GetDimensionSize(trial, "time", &num_timepoints);
	result &= data.GetDimensionSize(trial, "patients", &num_patients);
	result &= data.GetValues(trial, "patients", 0, num_patients, patient_ids);

	if (pk_type == PKMT_OneCompartment) {
		num_pk_params = 4;
		num_pk_pop_params = 4;
	} else if (pk_type == PKMT_TwoCompartment) {
		num_pk_params = 5;
		num_pk_pop_params = 4;
	} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
		num_pk_params = 7;
		num_pk_pop_params = 5;
	} else if (pk_type == PKMT_OneCompartmentTransit) {
		num_pk_params = 5;
		num_pk_pop_params = 4;
	} else if (pk_type == PKMT_TwoCompartmentTransit) {
		num_pk_params = 6;
		num_pk_pop_params = 4;
	} else {
		LOGERROR("Invalid PK model type");
		return false;
	}

	if (varset->GetNumVariables() != num_pk_params + num_pk_pop_params * (num_patients + 1) + 1) {
		LOGERROR("Incorrect number of variables in prior");
		return false;
	}

	time.setConstant(num_timepoints, std::numeric_limits<Real>::quiet_NaN());
	for (size_t i = 0; i < num_timepoints; i++) {
		result &= data.GetValue(trial, "time", i, time.data() + i);
	}

	Real minimum_dose = std::numeric_limits<Real>::max();
	observed_concentration.resize(num_patients);
	dose.resize(num_patients);
	dose_after_dose_change.resize(num_patients);
	dose_change_time.resize(num_patients);
	dosing_interval.resize(num_patients);
	intermittent.resize(num_patients, 0);
	skipped_days.resize(num_patients);
	simulate_until.resize(num_patients);
	for (size_t j = 0; j < num_patients; j++) {
		result &= data.GetValuesDim2(trial, drug + "_plasma_concentration", j, 0, num_timepoints, observed_concentration[j]);
		result &= data.GetValue(trial, drug + "_dose", j, &dose[j]);
		result &= data.GetValue(trial, drug + "_dose_after_dose_change", j, &dose_after_dose_change[j]);
		result &= data.GetValue(trial, drug + "_dose_change_time", j, &dose_change_time[j]);
		result &= data.GetValue(trial, drug + "_dosing_interval", j, &dosing_interval[j]);
		result &= data.GetValue(trial, drug + "_intermittent", j, &intermittent[j]);

		std::vector<unsigned int> interruptions;
		for (int i = 0; i < 29; i++) {
			unsigned int tmp;
			result &= data.GetValue(trial, "treatment_interruptions", j, (size_t)i, &tmp);
			if (tmp) {
				skipped_days[j].insert(i);
			}
		}

		// Patients marked as intermittent from day 2 are those where we do not know exactly when any interruptions occured
		// So we'll only simulate those for the very first day
		if (skipped_days[j].find(1) != skipped_days[j].end()) {
			for (ptrdiff_t i = 0; i < time.size(); i++) {
				if (time[i] >= 24.0) {
					simulate_until[j] = i;
					break;
				}
			}
		} else {
			simulate_until[j] = time.size();
		}

		for (ptrdiff_t i = 0; i < observed_concentration[j].size(); i++) {
			if (!std::isnan(observed_concentration[j](i))) {
				if (time[i] > 15 * 24) {
					// First measurement is too far into the start of treatment, can't trust the earlier estimates
					simulate_until[j] = 0.0;
				}
				break;
			}
		}

		if (!std::isnan(dose_after_dose_change[j])) {
			if (std::isnan(dose_change_time[j])) {
				LOGERROR("Patient %d has dose change, but time of dose change is not specified.", j);
				return false;
			}
			if (dose_change_time[j] / dosing_interval[j] != 0.0) {
				LOGERROR("Dose change time for patient  %d is not an exact multiple of the dosing interval.", j);
				return false;
			}
		}

		// Store minimum dose across patients for setting integration tolerances
		if (dose[j] < minimum_dose) {
			minimum_dose = dose[j];
		}
		if (!std::isnan(dose_after_dose_change[j]) && dose_after_dose_change[j] < minimum_dose) {
			minimum_dose = dose_after_dose_change[j];
		}
	}

	// Allocate structures for parallel evaluation
	for (size_t threadix = 0; threadix < sampling_threads; threadix++) {
		//solvers[threadix] = std::make_unique<ODESolverCVODE>();
		solvers[threadix] = std::make_unique<ODESolverDP5>();
		if (pk_type == PKMT_OneCompartment) {
			ODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_OneCompartment, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(2, NULL);
		} else if (pk_type == PKMT_TwoCompartment) {
			ODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartment, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			ODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartmentBiphasic, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_OneCompartmentTransit) {
			ODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_OneCompartmentTransit, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_TwoCompartmentTransit) {
			ODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartmentTransit, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(5, NULL);
		} else {
			LOGERROR("Invalid PK model type");
			return false;
		}
		solvers[threadix]->SetTolerance(5e-7f, minimum_dose * 5e-7f);
		parallel_data[threadix].simulated_concentrations.resize(num_patients);
		parallel_data[threadix].stored_trajectories.resize(num_patients);
	}

	prev_parameters.resize(parallel_data.size());
	prev_llh.resize(parallel_data.size());
	for (size_t i = 0; i < prev_parameters.size(); i++) {
		prev_parameters[i].resize(num_patients, VectorReal::Constant(8, 0));
		prev_llh[i].resize(num_patients, std::numeric_limits<Real>::quiet_NaN());
	}
	prev_ix.resize(num_patients, 0);

	return result;
}

bool LikelihoodPopPKTrajectory::PostInitialize()
{
	return true;
}

bool LikelihoodPopPKTrajectory::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	size_t sdix = varset->GetVariableIndex("standard_deviation");
	Real sd = varset->TransformVariable(sdix, values[sdix]);

	for (size_t j = 0; j < patient_ids.size(); j++) {
		ParallelData& pd = parallel_data[threadix];
		pd.dose = dose[j];
		pd.dosing_interval = dosing_interval[j];
		pd.dose_after_dose_change = dose_after_dose_change[j];
		pd.dose_change_time = dose_change_time[j];
		pd.intermittent = intermittent[j];
		pd.skipped_days = &skipped_days[j];
		
		pd.k_absorption  = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 0], values[0], values[num_pk_params + 0]));
		pd.k_excretion	 = bcm3::fastpow10(values[1]);
		pd.k_vod		 = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 2], values[3], values[num_pk_params + 2]));
		pd.k_elimination = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 1], values[2], values[num_pk_params + 1])) / pd.k_vod;
		if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic || pk_type == PKMT_TwoCompartmentTransit) {
			pd.k_intercompartmental = bcm3::fastpow10(values[4]) / pd.k_vod;
		}
		if (pk_type == PKMT_OneCompartmentTransit) {
			pd.k_transit = bcm3::fastpow10(values[4]);
		}
		if (pk_type == PKMT_TwoCompartmentTransit) {
			pd.k_transit = bcm3::fastpow10(values[5]);
		}
		if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			// No between-subject variation in biphasic time
			pd.k_biphasic_switch_time = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 3], values[5], values[num_pk_params + 3]));

			// The biphasic logic assumes that the switching time is always bigger than 0 and strictly less than the treatment interval
			pd.k_biphasic_switch_time = std::min(pd.k_biphasic_switch_time, pd.dosing_interval - 1e-2);
			pd.k_absorption2 = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 4], values[6], values[num_pk_params + 4]));
		}

		VectorReal parameters = VectorReal::Constant(8, 0);
		parameters(0) = pd.k_absorption;
		parameters(1) = pd.k_excretion;
		parameters(2) = pd.k_vod;
		parameters(3) = pd.k_elimination;
		if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic || pk_type == PKMT_TwoCompartmentTransit) {
			parameters(4) = pd.k_intercompartmental;
		}
		if (pk_type == PKMT_OneCompartmentTransit || pk_type == PKMT_TwoCompartmentTransit) {
			parameters(5) = pd.k_transit;
		} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			parameters(5) = pd.k_biphasic_switch_time;
			parameters(6) = pd.k_absorption2;
		}
		parameters(7) = sd;

		buffer_spinlock.lock();
		int which_exactly_equal = -1;
		for (size_t hi = 0; hi < prev_parameters.size(); hi++) {
			bool exactly_equal = true;
			for (int pi = 0; pi < parameters.size(); pi++) {
				if (parameters(pi) != prev_parameters[hi][j](pi)) {
					exactly_equal = false;
					break;
				}
			}
			if (exactly_equal) {
				which_exactly_equal = (int)hi;
				break;
			}
		}
		if (which_exactly_equal != -1) {
			logp += prev_llh[which_exactly_equal][j];
			buffer_spinlock.unlock();
			continue;
		} else {
			buffer_spinlock.unlock();
		}

		solvers[threadix]->SetUserData((void*)threadix);

		if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			pd.biphasic_switch = true;
			pd.current_dose_time = 0;
			solvers[threadix]->SetDiscontinuity(pd.k_biphasic_switch_time, boost::bind(&LikelihoodPopPKTrajectory::TreatmentCallbackBiphasic, this, _1, _2), (void*)threadix);
		} else {
			pd.current_dose_time = pd.dosing_interval;
			solvers[threadix]->SetDiscontinuity(pd.dosing_interval, boost::bind(&LikelihoodPopPKTrajectory::TreatmentCallback, this, _1, _2), (void*)threadix);
		}

		OdeReal initial_conditions[5];
		initial_conditions[0] = pd.dose;
		initial_conditions[1] = 0.0;
		initial_conditions[2] = 0.0;
		initial_conditions[3] = 0.0;
		initial_conditions[4] = 0.0;

		// The simulated trajectories are in mg; divide by volume of distribution in liters gives mg/l
		// Divide by mol weight gives mmol/l so 1e6 to go to nM (the units used in the datafile)
		Real MW;
		if (drug == "lapatinib") {
			MW = 581.06;
		} else if (drug == "dacomitinib") {
			MW = 469.95;
		} else if (drug == "afatinib") {
			MW = 485.94;
		} else if (drug == "trametinib") {
			MW = 615.404;
		} else if (drug == "mirdametinib") {
			MW = 482.19;
		} else if (drug == "selumetinib") {
			MW = 457.68;
		} else {
			LOGERROR("Unknown drug \"%s\"", drug.c_str());
			return false;
		}
		Real conversion = (1e6 / MW) / pd.k_vod;

		OdeVectorReal simulate_time = time.segment(0, simulate_until[j]);

		Real patient_logllh = 0.0;
		if (simulate_time.size() > 0) {
			if (!solvers[threadix]->Simulate(initial_conditions, simulate_time, pd.simulated_trajectories)) {
#if 0
				LOG("Patient %u failed", j);
				LOG(" k_absorption: %g", pd.k_absorption);
				LOG(" k_excretion: %g", pd.k_excretion);
				LOG(" k_vod: %g", pd.k_vod);
				LOG(" k_elimination: %g", pd.k_elimination);
#endif
				patient_logllh = -std::numeric_limits<Real>::infinity();
			} else {
				pd.simulated_concentrations[j] = pd.simulated_trajectories.row(1) * conversion;
				pd.stored_trajectories[j] = pd.simulated_trajectories;
				for (size_t i = 0; i < simulate_time.size(); i++) {
					Real x = conversion * pd.simulated_trajectories(1, i);
					Real y = observed_concentration[j](i);
					if (!std::isnan(y)) {
						patient_logllh += bcm3::LogPdfTnu4(x, y, sd);
					}
					if (std::isnan(x)) {
						//LOGERROR("NaN in trajectory for patient %zu, timepoint %zu", j, i);
						patient_logllh = -std::numeric_limits<Real>::infinity();
						break;
					}
				}
			}
		}

		logp += patient_logllh;

		buffer_spinlock.lock();
		prev_parameters[prev_ix[j]][j] = parameters;
		prev_llh[prev_ix[j]][j] = patient_logllh;
		prev_ix[j]++;
		if (prev_ix[j] >= prev_parameters.size()) {
			prev_ix[j] = 0;
		}
		buffer_spinlock.unlock();

		if (logp == -std::numeric_limits<Real>::infinity()) {
			break;
		}
	}

	return true;
}

bool LikelihoodPopPKTrajectory::CalculateDerivative_OneCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1];

	return true;
}

bool LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_intercompartmental * y[1] + pd.k_intercompartmental * y[2];
	dydt[2] = pd.k_intercompartmental * y[1] - pd.k_intercompartmental * y[2];

	return true;
}

bool LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartmentBiphasic(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	Real current_k_absorption;
	if (pd.biphasic_switch) {
		current_k_absorption = pd.k_absorption;
	} else {
		current_k_absorption = pd.k_absorption2;
	}

	dydt[0] = -(current_k_absorption + pd.k_excretion) * y[0];
	dydt[1] = current_k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_intercompartmental * y[1] + pd.k_intercompartmental * y[2];
	dydt[2] = pd.k_intercompartmental * y[1] - pd.k_intercompartmental * y[2];

	return true;
}

bool LikelihoodPopPKTrajectory::CalculateDerivative_OneCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_transit * y[2] - pd.k_elimination * y[1];
	dydt[2] = pd.k_absorption * y[0] - pd.k_transit * y[2];

	return true;
}

bool LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_transit * y[3] - pd.k_elimination * y[1] - pd.k_intercompartmental * y[1] + pd.k_intercompartmental * y[4];
	dydt[2] = pd.k_absorption * y[0] - pd.k_transit * y[2];
	dydt[3] = pd.k_transit * y[2] - pd.k_transit * y[3];
	dydt[4] = pd.k_intercompartmental * y[1] - pd.k_intercompartmental * y[4];

	return true;
}

Real LikelihoodPopPKTrajectory::TreatmentCallback(OdeReal t, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];
	pd.current_dose_time += pd.dosing_interval;
	bool give_treatment = true;
	int day = static_cast<int>(t / 24.0);
	if (pd.skipped_days->find(day) != pd.skipped_days->end()) {
		give_treatment = false;
	}
	if (pd.intermittent == 1) {
		Real time_in_week = t - 7.0 * 24.0 * floor(t / (7.0 * 24.0));
		if (time_in_week >= 5.0 * 24.0) {
			// No dose in day 6 and 7
			give_treatment = false;
		}
	} else if (pd.intermittent == 2) {
		Real time_in_treatment_course = t - 28.0 * 24.0 * floor(t / (28.0 * 24.0));
		if (time_in_treatment_course >= 21.0 * 24.0) {
			// No dose in day 22 through 28
			give_treatment = false;
		}
	} else if (pd.intermittent == 3) {
		Real time_in_week = t - 7.0 * 24.0 * floor(t / (7.0 * 24.0));
		if (time_in_week >= 4.0 * 24.0) {
			// No dose in day 5, 6 and 7
			give_treatment = false;
		}
	}
	Real dose = pd.dose;
	if (t >= pd.dose_change_time) {
		dose = pd.dose_after_dose_change;
	}
	if (give_treatment) {
		solvers[threadix]->set_y(0, solvers[threadix]->get_y(0) + dose);
	}
	return pd.current_dose_time;
}

Real LikelihoodPopPKTrajectory::TreatmentCallbackBiphasic(OdeReal t, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];
	if (pd.biphasic_switch) {
		pd.biphasic_switch = false;
		pd.current_dose_time += pd.dosing_interval;
		return pd.current_dose_time;
	} else {
		bool give_treatment = true;
		int day = static_cast<int>(t / 24.0);
		if (pd.skipped_days->find(day) != pd.skipped_days->end()) {
			give_treatment = false;
		}
		if (pd.intermittent == 1) {
			Real time_in_week = t - 7.0 * 24.0 * floor(t / (7.0 * 24.0));
			if (time_in_week >= 5.0 * 24.0) {
				// No dose in day 6 and 7
				give_treatment = false;
			}
		} else if (pd.intermittent == 2) {
			Real time_in_treatment_course = t - 28.0 * 24.0 * floor(t / (28.0 * 24.0));
			if (time_in_treatment_course >= 21.0 * 24.0) {
				// No dose in day 22 through 28
				give_treatment = false;
			}
		} else if (pd.intermittent == 3) {
			Real time_in_week = t - 7.0 * 24.0 * floor(t / (7.0 * 24.0));
			if (time_in_week >= 4.0 * 24.0) {
				// No dose in day 5, 6 and 7
				give_treatment = false;
			}
		}
		Real dose = pd.dose;
		if (t >= pd.dose_change_time) {
			dose = pd.dose_after_dose_change;
		}
		if (give_treatment) {
			solvers[threadix]->set_y(0, solvers[threadix]->get_y(0) + dose);
			pd.biphasic_switch = true;
			return pd.current_dose_time + pd.k_biphasic_switch_time;
		} else {
			pd.current_dose_time += pd.dosing_interval;
			return pd.current_dose_time;
		}
	}
}
