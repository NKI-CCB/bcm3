#include "Utils.h"
#include "LikelihoodPharmacokineticTrajectory.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"
#include "ODESolverCVODE.h"
#include "ODESolverDP5.h"

#include <boost/property_tree/xml_parser.hpp>
using namespace boost::placeholders;

typedef std::function<bool(OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;

LikelihoodPharmacokineticTrajectory::ParallelData::ParallelData()
	: dose(std::numeric_limits<Real>::quiet_NaN())
	, dosing_interval(std::numeric_limits<Real>::quiet_NaN())
	, dose_after_dose_change(std::numeric_limits<Real>::quiet_NaN())
	, dose_change_time(std::numeric_limits<Real>::quiet_NaN())
	, intermittent(0)
	, k_absorption(std::numeric_limits<Real>::quiet_NaN())
	, k_excretion(std::numeric_limits<Real>::quiet_NaN())
	, k_elimination(std::numeric_limits<Real>::quiet_NaN())
	, k_vod(std::numeric_limits<Real>::quiet_NaN())
	, k_periphery_fwd(std::numeric_limits<Real>::quiet_NaN())
	, k_periphery_bwd(std::numeric_limits<Real>::quiet_NaN())
	, k_transit(std::numeric_limits<Real>::quiet_NaN())
	, n_transit(std::numeric_limits<Real>::quiet_NaN())
	, k_biphasic_switch_time(std::numeric_limits<Real>::quiet_NaN())
	, k_absorption2(std::numeric_limits<Real>::quiet_NaN())
	, biphasic_switch(std::numeric_limits<Real>::quiet_NaN())
	, last_treatment(0.0)
	, current_dose_time(std::numeric_limits<Real>::quiet_NaN())
{
}

LikelihoodPharmacokineticTrajectory::LikelihoodPharmacokineticTrajectory(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, dose(std::numeric_limits<Real>::quiet_NaN())
	, dosing_interval(std::numeric_limits<Real>::quiet_NaN())
	, intermittent(false)
	, pk_type(PKMT_Undefined)
	, fixed_vod(std::numeric_limits<Real>::quiet_NaN())
	, fixed_periphery_fwd(std::numeric_limits<Real>::quiet_NaN())
	, fixed_periphery_bwd(std::numeric_limits<Real>::quiet_NaN())
{
	solvers.resize(sampling_threads);
	parallel_data.resize(sampling_threads);
}

LikelihoodPharmacokineticTrajectory::~LikelihoodPharmacokineticTrajectory()
{
}

void LikelihoodPharmacokineticTrajectory::SetDrugDosing(std::string drug, Real dose, Real interval, bool intermittent)
{
	this->drug = drug;
	this->dose = dose;
	this->dosing_interval = interval;
	this->intermittent = intermittent;
}

void LikelihoodPharmacokineticTrajectory::SetPatientID(std::string patient)
{
	this->patient_id = patient;
}

void LikelihoodPharmacokineticTrajectory::SetPKModelType(std::string type)
{
	if (type == "one") {
		pk_type = PKMT_OneCompartment;
	} else if (type == "two") {
		pk_type = PKMT_TwoCompartment;
	} else if (type == "one_biphasic_uptake") {
		pk_type = PKMT_TwoCompartmentBiphasicUptake;
	} else if (type == "two_biphasic_uptake") {
		pk_type = PKMT_TwoCompartmentBiphasicUptake;
	} else if (type == "one_transit") {
		pk_type = PKMT_OneCompartmentTransit;
	} else if (type == "two_transit") {
		pk_type = PKMT_TwoCompartmentTransit;
	} else {
		pk_type = PKMT_Undefined;
	}
}

bool LikelihoodPharmacokineticTrajectory::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;
	
	std::string pk_type_str;
	std::string trial;
	try {
		boost::property_tree::ptree& modelnode = likelihood_node.get_child("pk_model");
		drug = modelnode.get<std::string>("<xmlattr>.drug");
		pk_type_str = modelnode.get<std::string>("<xmlattr>.type");
		trial = modelnode.get<std::string>("<xmlattr>.trial");
		patient_id = modelnode.get<std::string>("<xmlattr>.patient", "");

		fixed_vod = modelnode.get<Real>("<xmlattr>.volume_of_distribution", std::numeric_limits<Real>::quiet_NaN());
		fixed_periphery_fwd = modelnode.get<Real>("<xmlattr>.k_periphery_fwd", std::numeric_limits<Real>::quiet_NaN());
		fixed_periphery_bwd = modelnode.get<Real>("<xmlattr>.k_periphery_bwd", std::numeric_limits<Real>::quiet_NaN());
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	std::string use_patient_str = vm["pk.patient"].as<std::string>();
	if (!use_patient_str.empty()) {
		patient_id = use_patient_str;
	}
	if (patient_id.empty()) {
		LOGERROR("Patient ID has not been specified in either the likelihood or as command-line option.");
		return false;
	}

	size_t fixed_var_count = 0;
	if (!std::isnan(fixed_vod)) fixed_var_count++;
	if (!std::isnan(fixed_periphery_fwd)) fixed_var_count++;
	if (!std::isnan(fixed_periphery_bwd)) fixed_var_count++;

	SetPKModelType(pk_type_str);
#if 0
	if (pk_type == PKMT_OneCompartment) {
		if (varset->GetNumVariables() != 6 - fixed_var_count) {
			LOGERROR("Incorrect number of variables in prior - should be 5 variables for a one-compartment model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartment) {
		if (varset->GetNumVariables() != 8 - fixed_var_count) {
			LOGERROR("Incorrect number of variables in prior - should be 7 variables for a two-compartment model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
		if (varset->GetNumVariables() != 10 - fixed_var_count) {
			LOGERROR("Incorrect number of variables in prior - should be 9 variables for a two-compartment biphasic model");
			return false;
		}
	} else if (pk_type == PKMT_OneCompartmentTransit) {
		if (varset->GetNumVariables() != 8 - fixed_var_count) {
			LOGERROR("Incorrect number of variables in prior - should be 6 variables for a one-compartment transit model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartmentTransit) {
		if (varset->GetNumVariables() != 10 - fixed_var_count) {
			LOGERROR("Incorrect number of variables in prior - should be 8 variables for a two-compartment transit model");
			return false;
		}
	} else {
		LOGERROR("Invalid PK model type");
		return false;
	}
#endif

	bool result = true;

	bcm3::NetCDFDataFile data;
	if (!data.Open("pkdata.nc", false)) {
		return false;
	}

	size_t num_timepoints;
	result &= data.GetDimensionSize(trial, "time", &num_timepoints);

	// Find patient index
	size_t patient_ix;
	if (!data.GetDimensionIx(trial, "patients", patient_id, &patient_ix)) {
		LOGERROR("Cannot find patient \"%s\" in data file", patient_id.c_str());
		return false;
	}

	time.resize(num_timepoints);
	observed_concentration.resize(num_timepoints);
	for (size_t i = 0; i < num_timepoints; i++) {
		result &= data.GetValue(trial, "time", i, time.data() + i);
		result &= data.GetValue(trial, drug + "_plasma_concentration", patient_ix, i, observed_concentration.data() + i);
	}

	result &= data.GetValue(trial, drug + "_dose", patient_ix, &dose);
	result &= data.GetValue(trial, drug + "_dosing_interval", patient_ix, &dosing_interval);
	result &= data.GetValue(trial, drug + "_dose_after_dose_change", patient_ix, &dose_after_dose_change);
	result &= data.GetValue(trial, drug + "_dose_change_time", patient_ix, &dose_change_time);
	unsigned int intermittint_int;
	result &= data.GetValue(trial, drug + "_intermittent", patient_ix, &intermittint_int);
	intermittent = intermittint_int ? true : false;

	std::vector<unsigned int> interruptions;
	for (int i = 0; i < 29; i++) {
		unsigned int tmp;
		result &= data.GetValue(trial, "treatment_interruptions", patient_ix, (size_t)i, &tmp);
		if (tmp) {
			skipped_days.insert(i);
		}
	}

	// Allocate structures for parallel evaluation
	for (size_t threadix = 0; threadix < sampling_threads; threadix++) {
		//solvers[threadix] = std::make_unique<ODESolverDP5>();
		solvers[threadix] = std::make_unique<ODESolverCVODE>();
		if (pk_type == PKMT_OneCompartment) {
			solvers[threadix]->SetDerivativeFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartment, this, _1, _2, _3, _4));
			solvers[threadix]->SetJacobianFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateJacobian_OneCompartment, this, _1, _2, _3, _4, _5));
			solvers[threadix]->Initialize(2, NULL);
		} else if (pk_type == PKMT_TwoCompartment) {
			solvers[threadix]->SetDerivativeFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartment, this, _1, _2, _3, _4));
			solvers[threadix]->SetJacobianFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateJacobian_TwoCompartment, this, _1, _2, _3, _4, _5));
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_OneCompartmentBiphasicUptake) {
			solvers[threadix]->SetJacobianFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateJacobian_OneCompartmentBiphasicUptake, this, _1, _2, _3, _4, _5));
			solvers[threadix]->SetDerivativeFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartmentBiphasicUptake, this, _1, _2, _3, _4));
			solvers[threadix]->Initialize(2, NULL);
		} else if (pk_type == PKMT_TwoCompartmentBiphasicUptake) {
			solvers[threadix]->SetJacobianFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateJacobian_TwoCompartmentBiphasicUptake, this, _1, _2, _3, _4, _5));
			solvers[threadix]->SetDerivativeFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentBiphasicUptake, this, _1, _2, _3, _4));
			solvers[threadix]->Initialize(3, NULL);
		} else if (pk_type == PKMT_OneCompartmentTransit) {
			solvers[threadix]->SetDerivativeFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartmentTransit, this, _1, _2, _3, _4));
			solvers[threadix]->SetJacobianFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateJacobian_OneCompartmentTransit, this, _1, _2, _3, _4, _5));
			solvers[threadix]->Initialize(2, NULL);
		} else if (pk_type == PKMT_TwoCompartmentTransit) {
			solvers[threadix]->SetDerivativeFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentTransit, this, _1, _2, _3, _4));
			solvers[threadix]->SetJacobianFunction(boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateJacobian_TwoCompartmentTransit, this, _1, _2, _3, _4, _5));
			solvers[threadix]->Initialize(3, NULL);
		} else {
			LOGERROR("Invalid PK model type");
			return false;
		}
		solvers[threadix]->SetTolerance(1e-6f, dose * 1e-6f);
	}

	return result;
}

bool LikelihoodPharmacokineticTrajectory::PostInitialize()
{
	return true;
}

bool LikelihoodPharmacokineticTrajectory::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	size_t sdix = varset->GetVariableIndex("standard_deviation");
	Real sd = varset->TransformVariable(sdix, values[sdix]);
	Real sd2 = varset->TransformVariable(sdix+1, values[sdix+1]);

	size_t absorption_ix = 0;
	size_t excretion_ix = 1;
	size_t elimination_ix = 2;
	size_t vod_ix = 3;
	size_t periphery_fwd_ix = 4;
	size_t periphery_bwd_ix = 5;

	ParallelData& pd = parallel_data[threadix];
	pd.dose				= dose;
	pd.intermittent		= intermittent;
	pd.dosing_interval	= dosing_interval;
	pd.dose_after_dose_change = dose_after_dose_change;
	pd.dose_change_time = dose_change_time;
	pd.k_absorption		= varset->TransformVariable(absorption_ix, values[absorption_ix]);
	pd.k_excretion		= varset->TransformVariable(excretion_ix, values[excretion_ix]);
	pd.k_vod			= std::isnan(fixed_vod) ? varset->TransformVariable(vod_ix, values[vod_ix]) : fixed_vod;
	pd.k_elimination	= varset->TransformVariable(elimination_ix, values[elimination_ix]) / pd.k_vod;
	if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiphasicUptake || pk_type == PKMT_TwoCompartmentTransit) {
		if (std::isnan(fixed_periphery_fwd)) {
			pd.k_periphery_fwd = varset->TransformVariable(periphery_fwd_ix, values[periphery_fwd_ix]);
			pd.k_periphery_bwd = varset->TransformVariable(periphery_bwd_ix, values[periphery_bwd_ix]);
		} else {
			pd.k_periphery_fwd = fixed_periphery_fwd;
			pd.k_periphery_bwd = fixed_periphery_bwd;
		}
	}
	if (pk_type == PKMT_OneCompartmentTransit || pk_type == PKMT_TwoCompartmentTransit) {
		size_t n_transit_ix = varset->GetVariableIndex("n_transit");
		pd.n_transit = varset->TransformVariable(n_transit_ix, values[n_transit_ix]);
		size_t transit_time_ix = varset->GetVariableIndex("mean_transit_time");
		pd.k_transit = (pd.n_transit + 1) / varset->TransformVariable(transit_time_ix, values[transit_time_ix]);
	}
	if (pk_type == PKMT_OneCompartmentBiphasicUptake || pk_type == PKMT_TwoCompartmentBiphasicUptake) {
		// The biphasic logic assumes that the switching time is always bigger than 0 and strictly less than the treatment interval
		pd.k_biphasic_switch_time = varset->TransformVariable(6, values[6]);
		pd.k_absorption2 = varset->TransformVariable(7, values[7]);
		pd.biphasic_switch = true;
		pd.current_dose_time = 0;
		solvers[threadix]->SetDiscontinuity(pd.k_biphasic_switch_time, boost::bind(&LikelihoodPharmacokineticTrajectory::TreatmentCallbackBiphasic, this, _1, _2), (void*)threadix);
	} else {
		pd.current_dose_time = dosing_interval;
		solvers[threadix]->SetDiscontinuity(dosing_interval, boost::bind(&LikelihoodPharmacokineticTrajectory::TreatmentCallback, this, _1, _2), (void*)threadix);
	}
	pd.skipped_days = &skipped_days;
	pd.last_treatment = 0.0;
	solvers[threadix]->SetUserData((void*)threadix);

	OdeReal initial_conditions[3];
	if (pk_type == PKMT_OneCompartmentTransit || pk_type == PKMT_TwoCompartmentTransit) {
		initial_conditions[0] = 0.0;
	} else {
		initial_conditions[0] = dose;
	}
	initial_conditions[1] = 0.0;
	initial_conditions[2] = 0.0;

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

	if (!solvers[threadix]->Simulate(initial_conditions, time, pd.simulated_trajectories)) {
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	} else {
		pd.simulated_concentrations = pd.simulated_trajectories.row(1) * conversion;

		for (size_t i = 0; i < time.size(); i++) {
			Real x = pd.simulated_concentrations(i);
			Real y = observed_concentration[i];
			if (!std::isnan(y)) {
				logp += bcm3::LogPdfTnu4(x, y, sd + sd2 * std::max(x, 0.0));
			}
		}
	}

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateJacobian_OneCompartment(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	jac(0, 0) = -(pd.k_absorption + pd.k_excretion);
	jac(1, 0) = pd.k_absorption;
	jac(1, 1) = -pd.k_elimination;

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_periphery_fwd * y[1] + pd.k_periphery_bwd * y[2];
	dydt[2] = pd.k_periphery_fwd * y[1] - pd.k_periphery_bwd * y[2];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateJacobian_TwoCompartment(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	jac(0, 0) = -(pd.k_absorption + pd.k_excretion);
	jac(1, 0) = pd.k_absorption;
	jac(1, 1) = -(pd.k_elimination + pd.k_periphery_fwd);
	jac(1, 2) = pd.k_periphery_bwd;
	jac(2, 1) = pd.k_periphery_fwd;
	jac(2, 2) = -pd.k_periphery_bwd;

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
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
	dydt[1] = current_k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_periphery_fwd * y[1];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateJacobian_OneCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	Real current_k_absorption;
	if (pd.biphasic_switch) {
		current_k_absorption = pd.k_absorption;
	} else {
		current_k_absorption = pd.k_absorption2;
	}

	jac(0, 0) = -(current_k_absorption + pd.k_excretion);
	jac(1, 0) = current_k_absorption;
	jac(1, 1) = -pd.k_elimination;

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
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
	dydt[1] = current_k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_periphery_fwd * y[1] + pd.k_periphery_bwd * y[2];
	dydt[2] = pd.k_periphery_fwd * y[1] - pd.k_periphery_bwd * y[2];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateJacobian_TwoCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	Real current_k_absorption;
	if (pd.biphasic_switch) {
		current_k_absorption = pd.k_absorption;
	} else {
		current_k_absorption = pd.k_absorption2;
	}

	jac(0, 0) = -(current_k_absorption + pd.k_excretion);
	jac(1, 0) = current_k_absorption;
	jac(1, 1) = -(pd.k_elimination + pd.k_periphery_fwd);
	jac(1, 2) = pd.k_periphery_bwd;
	jac(2, 1) = pd.k_periphery_fwd;
	jac(2, 2) = -pd.k_periphery_bwd;

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	Real dose = pd.dose;
	if (t >= pd.dose_change_time) {
		dose = pd.dose_after_dose_change;
	}

	Real t_since_treatment = t - pd.last_treatment;
	Real log_n_transit_factorial = 0.9189385332046727 + (pd.n_transit + 0.5) * log(pd.n_transit) - pd.n_transit + log(1 + 1 / (12.0 * pd.n_transit));
	Real transit = exp((pd.n_transit * log(pd.k_transit * t_since_treatment) - pd.k_transit * t_since_treatment) - log_n_transit_factorial);
	transit = pd.k_transit * transit * dose;

	dydt[0] = transit - (pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateJacobian_OneCompartmentTransit(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	jac(0, 0) = -(pd.k_absorption + pd.k_excretion);
	jac(1, 0) = pd.k_absorption;
	jac(1, 1) = -pd.k_elimination;

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	Real dose = pd.dose;
	if (t >= pd.dose_change_time) {
		dose = pd.dose_after_dose_change;
	}

	Real t_since_treatment = t - pd.last_treatment;
	Real log_n_transit_factorial = 0.9189385332046727 + (pd.n_transit + 0.5) * log(pd.n_transit) - pd.n_transit + log(1 + 1 / (12.0 * pd.n_transit));
	Real transit = exp((pd.n_transit * log(pd.k_transit * t_since_treatment) - pd.k_transit * t_since_treatment) - log_n_transit_factorial);
	transit = pd.k_transit * transit * dose;

	dydt[0] = transit - (pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_periphery_fwd * y[1] + pd.k_periphery_bwd * y[2];
	dydt[2] = pd.k_periphery_fwd * y[1] - pd.k_periphery_bwd * y[2];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateJacobian_TwoCompartmentTransit(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	jac(0, 0) = -(pd.k_absorption + pd.k_excretion);
	jac(1, 0) = pd.k_absorption;
	jac(1, 1) = -(pd.k_elimination + pd.k_periphery_fwd);
	jac(1, 2) = pd.k_periphery_bwd;
	jac(2, 1) = pd.k_periphery_fwd;
	jac(2, 2) = -pd.k_periphery_bwd;

	return true;
}

inline bool LikelihoodPharmacokineticTrajectory::CheckGiveTreatment(OdeReal t, ParallelData& pd)
{
	bool give_treatment = true;
	int day = static_cast<int>(floor(t / 24.0));
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
	return give_treatment;
}

Real LikelihoodPharmacokineticTrajectory::TreatmentCallback(OdeReal t, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];
	pd.current_dose_time += pd.dosing_interval;
	if (CheckGiveTreatment(t, pd)) {
		Real dose = pd.dose;
		if (t >= pd.dose_change_time) {
			dose = pd.dose_after_dose_change;
		}
		if (pk_type == PKMT_OneCompartmentTransit || pk_type == PKMT_TwoCompartmentTransit) {
			pd.last_treatment = t;
		} else {
			solvers[threadix]->set_y(0, solvers[threadix]->get_y(0) + dose);
		}
	}
	return pd.current_dose_time;
}

Real LikelihoodPharmacokineticTrajectory::TreatmentCallbackBiphasic(OdeReal t, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];
	if (pd.biphasic_switch) {
		pd.biphasic_switch = false;
		pd.current_dose_time += pd.dosing_interval;
		return pd.current_dose_time;
	} else {
		if (CheckGiveTreatment(t, pd)) {
			Real dose = pd.dose;
			if (t >= pd.dose_change_time) {
				dose = pd.dose_after_dose_change;
			}
			if (pk_type == PKMT_OneCompartmentTransit || pk_type == PKMT_TwoCompartmentTransit) {
				pd.last_treatment = t;
			} else {
				solvers[threadix]->set_y(0, solvers[threadix]->get_y(0) + dose);
			}
			pd.biphasic_switch = true;
			return pd.current_dose_time + pd.k_biphasic_switch_time;
		} else {
			pd.current_dose_time += pd.dosing_interval;
			return pd.current_dose_time;
		}
	}
}

void LikelihoodPharmacokineticTrajectory::AddOptionsDescription(boost::program_options::options_description& pod)
{
	pod.add_options()
		("pk.patient", boost::program_options::value<std::string>()->default_value(""), "Fit the data from this specific patient.")
	;
}
