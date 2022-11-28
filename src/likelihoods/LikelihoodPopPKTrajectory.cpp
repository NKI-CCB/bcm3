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
	, dose_cycle2(std::numeric_limits<Real>::quiet_NaN())
	, dosing_interval(std::numeric_limits<Real>::quiet_NaN())
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

	num_patients = 3;

	result &= data.GetValues(trial, "patients", 0, num_patients, patient_ids);

	if (pk_type == PKMT_OneCompartment) {
		num_pk_params = 4;
		num_pk_pop_params = 4;
		if (varset->GetNumVariables() != num_pk_params + num_pk_pop_params * (num_patients + 1) + 1) {
			LOGERROR("Incorrect number of variables in prior - should be 4 variables per patient for a one-compartment model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartment) {
		num_pk_params = 6;
		num_pk_pop_params = 6;
		if (varset->GetNumVariables() != num_pk_params + num_pk_pop_params * (num_patients + 1) + 1) {
			LOGERROR("Incorrect number of variables in prior - should be 6 variables for a two-compartment model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
		num_pk_params = 7;
		num_pk_pop_params = 5;
		if (varset->GetNumVariables() != num_pk_params + num_pk_pop_params * (num_patients + 1) + 1) {
			LOGERROR("Incorrect number of variables in prior - should be 7 mean variables and 5 per-patient variables for a two-compartment biphasic model");
			return false;
		}
	} else if (pk_type == PKMT_OneCompartmentTransit) {
		num_pk_params = 5;
		num_pk_pop_params = 5;
		if (varset->GetNumVariables() != num_pk_params + num_pk_pop_params * (num_patients + 1) + 1) {
			LOGERROR("Incorrect number of variables in prior - should be 5 variables for a one-compartment transit model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartmentTransit) {
		num_pk_params = 7;
		num_pk_pop_params = 7;
		if (varset->GetNumVariables() != num_pk_params + num_pk_pop_params * (num_patients + 1) + 1) {
			LOGERROR("Incorrect number of variables in prior - should be 7 variables for a two-compartment transit model");
			return false;
		}
	} else {
		LOGERROR("Invalid PK model type");
		return false;
	}

	result &= data.GetValues(trial, "time", 0, num_timepoints, time);

	Real minimum_dose = std::numeric_limits<Real>::max();
	observed_concentration.resize(num_patients);
	dose.resize(num_patients);
	dose_cycle2.resize(num_patients);
	dosing_interval.resize(num_patients);
	intermittent.resize(num_patients, 0);
	skipped_days.resize(num_patients);
	for (size_t j = 0; j < num_patients; j++) {
		observed_concentration[j].resize(num_timepoints);
		for (size_t i = 0; i < num_timepoints; i++) {
			result &= data.GetValue(trial, drug + "_plasma_concentration", j, i, observed_concentration[j].data() + i);
		}

		result &= data.GetValue(trial, drug + "_dose", j, &dose[j]);
		result &= data.GetValue(trial, drug + "_dose_cycle2", j, &dose_cycle2[j]);
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

		if (std::isnan(dose_cycle2[j])) {
			dose_cycle2[j] = dose[j];
		}

		if (dose[j] < minimum_dose) {
			minimum_dose = dose[j];
		}
		if (dose_cycle2[j] < minimum_dose) {
			minimum_dose = dose_cycle2[j];
		}
	}

	// Allocate structures for parallel evaluation
	for (size_t threadix = 0; threadix < sampling_threads; threadix++) {
		//solvers[threadix] = std::make_unique<CVODESolver>();
		solvers[threadix] = std::make_unique<ODESolverDP5>();
		if (pk_type == PKMT_OneCompartment) {
			CVODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_OneCompartment, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(2, NULL);
		} else if (pk_type == PKMT_TwoCompartment) {
			CVODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartment, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			CVODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartmentBiphasic, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_OneCompartmentTransit) {
			CVODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_OneCompartmentTransit, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_TwoCompartmentTransit) {
			CVODESolver::TDeriviativeFunction derivative = boost::bind(&LikelihoodPopPKTrajectory::CalculateDerivative_TwoCompartmentTransit, this, _1, _2, _3, _4);
			solvers[threadix]->SetDerivativeFunction(derivative);
			solvers[threadix]->Initialize(5, NULL);
		} else {
			LOGERROR("Invalid PK model type");
			return false;
		}
		solvers[threadix]->SetTolerance(1e-8, minimum_dose * 1e-8);
		parallel_data[threadix].simulated_concentrations.resize(num_patients);
		parallel_data[threadix].stored_trajectories.resize(num_patients);
	}

	prev_parameters.resize(parallel_data.size());
	prev_llh.resize(parallel_data.size());
	for (size_t i = 0; i < prev_parameters.size(); i++) {
		prev_parameters[i].resize(num_patients, VectorReal::Constant(10, 0));
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
		pd.dose_cycle2 = dose_cycle2[j];
		pd.intermittent = intermittent[j];
		pd.dosing_interval = dosing_interval[j];
		pd.skipped_days = &skipped_days[j];
		
#if 0
		pd.k_absorption  =  varset->TransformVariable(0, values[0]);
		pd.k_excretion   =  varset->TransformVariable(1, values[1]);
		pd.k_vod		 =  varset->TransformVariable(3, values[3]);
		pd.k_elimination = (varset->TransformVariable(2, values[2])) / pd.k_vod;
		if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic) {
			pd.k_intercompartmental = varset->TransformVariable(4, values[4]);
			pd.k_intercompartmental = varset->TransformVariable(5, values[5]);
		}
		if (pk_type == PKMT_OneCompartmentTransit) {
			pd.k_transit = varset->TransformVariable(4, values[4]);
		}
		if (pk_type == PKMT_TwoCompartmentTransit) {
			pd.k_transit = varset->TransformVariable(6, values[6]);
		}
#elif 0
		pd.k_absorption  =  varset->TransformVariable(0, values[0]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 0], 0, exp(values[num_pk_params + 0])));
		pd.k_excretion   =  varset->TransformVariable(1, values[1]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 1], 0, exp(values[num_pk_params + 1])));
		pd.k_vod		 =  varset->TransformVariable(3, values[3]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 3], 0, exp(values[num_pk_params + 3])));
		pd.k_elimination = (varset->TransformVariable(2, values[2]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 2], 0, exp(values[num_pk_params + 2])))) / pd.k_vod;
		if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic || pk_type == PKMT_TwoCompartmentTransit) {
			pd.k_intercompartmental = varset->TransformVariable(4, values[4]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 4], 0, exp(values[num_pk_params + 4])));
			pd.k_intercompartmental = varset->TransformVariable(5, values[5]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 5], 0, exp(values[num_pk_params + 5])));
		}
		if (pk_type == PKMT_OneCompartmentTransit) {
			pd.k_transit = varset->TransformVariable(4, values[4]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 4], 0, exp(values[num_pk_params + 4])));
		}
		if (pk_type == PKMT_TwoCompartmentTransit) {
			pd.k_transit = varset->TransformVariable(6, values[6]) * exp(bcm3::QuantileNormal(values[num_pk_params * (j + 2) + 6], 0, exp(values[num_pk_params + 6])));
		}
#else
		pd.k_absorption  = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 0], values[0], values[num_pk_params + 0]));
		//pd.k_excretion   = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 1], values[1], values[num_pk_params + 1]));
		pd.k_excretion = bcm3::fastpow10(values[1]);
		pd.k_vod		 = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 2], values[3], values[num_pk_params + 2]));
		pd.k_elimination = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 1], values[2], values[num_pk_params + 1])) / pd.k_vod;
		if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic || pk_type == PKMT_TwoCompartmentTransit) {
			//pd.k_intercompartmental = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 4], values[4], values[num_pk_params + 4]));
			//pd.k_intercompartmental = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 5], values[5], values[num_pk_params + 5]));
			pd.k_intercompartmental = bcm3::fastpow10(values[4]) / pd.k_vod;
		}
		if (pk_type == PKMT_OneCompartmentTransit) {
			pd.k_transit = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 4], values[4], values[num_pk_params + 4]));
		}
		if (pk_type == PKMT_TwoCompartmentTransit) {
			pd.k_transit = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 6], values[6], values[num_pk_params + 6]));
		}
		if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			// No between-subject variation in biphasic time
			pd.k_biphasic_switch_time = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 3], values[5], values[num_pk_params + 3]));

			// The biphasic logic assumes that the switching time is always bigger than 0 and strictly less than the treatment interval
			pd.k_biphasic_switch_time = std::min(pd.k_biphasic_switch_time, pd.dosing_interval - 1e-2);
			pd.k_absorption2 = bcm3::fastpow10(bcm3::QuantileNormal(values[num_pk_params + num_pk_pop_params * (j + 1) + 4], values[6], values[num_pk_params + 4]));
		}
#endif

		VectorReal parameters = VectorReal::Constant(10, 0);
		parameters(0) = pd.k_absorption;
		parameters(1) = pd.k_excretion;
		parameters(2) = pd.k_vod;
		parameters(3) = pd.k_elimination;
		if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic || pk_type == PKMT_TwoCompartmentTransit) {
			parameters(4) = pd.k_intercompartmental;
		}
		if (pk_type == PKMT_OneCompartmentTransit) {
			parameters(5) = pd.k_transit;
		}
		if (pk_type == PKMT_TwoCompartmentTransit) {
			parameters(5) = pd.k_transit;
		}
		if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			parameters(6) = pd.k_biphasic_switch_time;
			parameters(7) = pd.k_absorption2;
		}
		parameters(8) = sd;

		buffer_spinlock.lock();
		int which_exactly_equal = -1;
		for (size_t hi = 0; hi < prev_parameters.size(); hi++) {
			bool exactly_equal = true;
			for (int pi = 0; pi < 10; pi++) {
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

		Real initial_conditions[5];
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

		Real patient_logllh = 0.0;
		if (!solvers[threadix]->Simulate(initial_conditions, time, pd.simulated_trajectories)) {
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
			for (size_t i = 0; i < time.size(); i++) {
				Real x = conversion * pd.simulated_trajectories(1, i);
				Real y = observed_concentration[j](i);
				if (y == y) {
					patient_logllh += bcm3::LogPdfTnu3(x, y, sd);
				}
				if (std::isnan(x)) {
					//LOGERROR("NaN in trajectory for patient %zu, timepoint %zu", j, i);
					patient_logllh = -std::numeric_limits<Real>::infinity();
					break;
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
	if (day >= 28) {
		dose = pd.dose_cycle2;
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
		if (day >= 28) {
			dose = pd.dose_cycle2;
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
