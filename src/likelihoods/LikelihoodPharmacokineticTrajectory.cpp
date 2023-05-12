#include "Utils.h"
#include "LikelihoodPharmacokineticTrajectory.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include <boost/property_tree/xml_parser.hpp>
using namespace boost::placeholders;

typedef std::function<bool(OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;


LikelihoodPharmacokineticTrajectory::LikelihoodPharmacokineticTrajectory(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, dose(std::numeric_limits<Real>::quiet_NaN())
	, dosing_interval(std::numeric_limits<Real>::quiet_NaN())
	, intermittent(false)
	, pk_type(PKMT_Undefined)
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
	} else if (type == "two_biphasic") {
		pk_type = PKMT_TwoCompartmentBiPhasic;
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

	if (pk_type == PKMT_OneCompartment) {
		if (varset->GetNumVariables() != 5) {
			LOGERROR("Incorrect number of variables in prior - should be 5 variables for a one-compartment model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartment) {
		if (varset->GetNumVariables() != 7) {
			LOGERROR("Incorrect number of variables in prior - should be 7 variables for a two-compartment model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
		if (varset->GetNumVariables() != 9) {
			LOGERROR("Incorrect number of variables in prior - should be 9 variables for a two-compartment biphasic model");
			return false;
		}
	} else if (pk_type == PKMT_OneCompartmentTransit) {
		if (varset->GetNumVariables() != 6) {
			LOGERROR("Incorrect number of variables in prior - should be 6 variables for a one-compartment transit model");
			return false;
		}
	} else if (pk_type == PKMT_TwoCompartmentTransit) {
		if (varset->GetNumVariables() != 8) {
			LOGERROR("Incorrect number of variables in prior - should be 8 variables for a two-compartment transit model");
			return false;
		}
	} else {
		LOGERROR("Invalid PK model type");
		return false;
	}

	bool result = true;

	bcm3::NetCDFDataFile data;
	if (!data.Open("pkdata.nc", false)) {
		return false;
	}

	size_t num_timepoints;
	result &= data.GetDimensionSize("LTK", "time", &num_timepoints);

	// Find patient index
	size_t patient_ix;
	data.GetDimensionIx("LTK", "patients", patient_id, &patient_ix);

	time.resize(num_timepoints);
	observed_concentration.resize(num_timepoints);
	for (size_t i = 0; i < num_timepoints; i++) {
		result &= data.GetValue("LTK", "time", i, time.data() + i);
		result &= data.GetValue("LTK", drug, patient_ix, i, observed_concentration.data() + i);
	}

	// Allocate structures for parallel evaluation
	for (size_t threadix = 0; threadix < sampling_threads; threadix++) {
		if (pk_type == PKMT_OneCompartment) {
			ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartment, this, _1, _2, _3, _4);
			solvers[threadix].SetDerivativeFunction(derivative);
			solvers[threadix].Initialize(2, NULL);
		} else if (pk_type == PKMT_TwoCompartment) {
			ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartment, this, _1, _2, _3, _4);
			solvers[threadix].SetDerivativeFunction(derivative);
			solvers[threadix].Initialize(3, NULL);
		} else if (pk_type == PKMT_TwoCompartmentBiPhasic) {
			ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentBiphasic, this, _1, _2, _3, _4);
			solvers[threadix].SetDerivativeFunction(derivative);
			solvers[threadix].Initialize(3, NULL);
		} else if (pk_type == PKMT_OneCompartmentTransit) {
			ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartmentTransit, this, _1, _2, _3, _4);
			solvers[threadix].SetDerivativeFunction(derivative);
			solvers[threadix].Initialize(3, NULL);
		} else if (pk_type == PKMT_TwoCompartmentTransit) {
			ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentTransit, this, _1, _2, _3, _4);
			solvers[threadix].SetDerivativeFunction(derivative);
			solvers[threadix].Initialize(5, NULL);
		} else {
			LOGERROR("Invalid PK model type");
			return false;
		}
		solvers[threadix].SetTolerance(1e-6, 1e-6);
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

	ParallelData& pd = parallel_data[threadix];
	pd.k_absorption		= varset->TransformVariable(0, values[0]);
	pd.k_excretion		= varset->TransformVariable(1, values[1]);
	pd.k_elimination	= varset->TransformVariable(2, values[2]) / varset->TransformVariable(3, values[3]);
	if (pk_type == PKMT_TwoCompartment || pk_type == PKMT_TwoCompartmentBiPhasic) {
		pd.k_periphery_fwd = varset->TransformVariable(5, values[5]);
		pd.k_periphery_bwd = varset->TransformVariable(6, values[6]);
	}
	if (pk_type == PKMT_OneCompartmentTransit) {
		pd.k_transit = varset->TransformVariable(5, values[5]);
	}
	if (pk_type == PKMT_TwoCompartmentTransit) {
		pd.k_transit = varset->TransformVariable(7, values[7]);
	}
	solvers[threadix].SetUserData((void*)threadix);
	if (pk_type == PKMT_TwoCompartmentBiPhasic) {
		// The biphasic logic assumes that the switching time is always bigger than 0 and strictly less than the treatment interval
		pd.k_biphasic_switch_time = varset->TransformVariable(7, values[7]);
		pd.k_absorption2 = varset->TransformVariable(8, values[8]);
		pd.biphasic_switch = true;
		pd.current_dose_time = 0;
		solvers[threadix].SetDiscontinuity(pd.k_biphasic_switch_time, boost::bind(&LikelihoodPharmacokineticTrajectory::TreatmentCallbackBiphasic, this, _1, _2), (void*)threadix);
	} else {
		pd.current_dose_time = dosing_interval;
		solvers[threadix].SetDiscontinuity(dosing_interval, boost::bind(&LikelihoodPharmacokineticTrajectory::TreatmentCallback, this, _1, _2), (void*)threadix);
	}

	OdeReal initial_conditions[5];
	initial_conditions[0] = dose;
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
	} else if (drug == "trametinib") {
		MW = 615.404;
	} else {
		LOGERROR("Unknown drug \"%s\"", drug.c_str());
		return false;
	}
	Real conversion = (1e6 / MW) / varset->TransformVariable(3, values[3]);

	if (!solvers[threadix].Simulate(initial_conditions, time, pd.simulated_trajectories)) {
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	} else {
		size_t count = 0;
		for (size_t i = 0; i < time.size(); i++) {
			Real x = conversion * pd.simulated_trajectories(1, i);
			Real y = observed_concentration[i];
			if (y == y) {
				logp += bcm3::LogPdfTnu3(x, y, values[4]);
				count++;
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

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_absorption * y[0] - pd.k_elimination * y[1] - pd.k_periphery_fwd * y[1] + pd.k_periphery_bwd * y[2];
	dydt[2] = pd.k_periphery_fwd * y[1] - pd.k_periphery_bwd * y[2];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentBiphasic(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
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

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_OneCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_transit * y[2] - pd.k_elimination * y[1];
	dydt[2] = pd.k_absorption * y[0] - pd.k_transit * y[2];

	return true;
}

bool LikelihoodPharmacokineticTrajectory::CalculateDerivative_TwoCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];

	dydt[0] = -(pd.k_absorption + pd.k_excretion) * y[0];
	dydt[1] = pd.k_transit * y[3] - pd.k_elimination * y[1] - pd.k_periphery_fwd * y[1] + pd.k_periphery_bwd * y[4];
	dydt[2] = pd.k_absorption * y[0] - pd.k_transit * y[2];
	dydt[3] = pd.k_transit * y[2] - pd.k_transit * y[3];
	dydt[4] = pd.k_periphery_fwd * y[1] - pd.k_periphery_bwd * y[4];

	return true;
}

Real LikelihoodPharmacokineticTrajectory::TreatmentCallback(OdeReal t, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];
	pd.current_dose_time += dosing_interval;
	if (intermittent) {
		Real time_in_week = t - floor(t / 7.0 * 24.0);
		if (time_in_week < 5.0 * 24.0) {
			solvers[threadix].set_y(0, solvers[threadix].get_y(0) + dose);
		} else {
			// No dose in day 6 and 7
		}
	} else {
		solvers[threadix].set_y(0, solvers[threadix].get_y(0) + dose);
	}
	return pd.current_dose_time;
}

Real LikelihoodPharmacokineticTrajectory::TreatmentCallbackBiphasic(OdeReal t, void* user)
{
	size_t threadix = (size_t)user;
	ParallelData& pd = parallel_data[threadix];
	if (pd.biphasic_switch) {
		pd.biphasic_switch = false;
		pd.current_dose_time += dosing_interval;
		return pd.current_dose_time;
	} else {
		pd.biphasic_switch = true;

		if (intermittent) {
			Real time_in_week = t - floor(t / 7.0 * 24.0);
			if (time_in_week < 5.0 * 24.0) {
				solvers[threadix].set_y(0, solvers[threadix].get_y(0) + dose);
			} else {
				// No dose in day 6 and 7
			}
		} else {
			solvers[threadix].set_y(0, solvers[threadix].get_y(0) + dose);
		}
		return pd.current_dose_time + pd.k_biphasic_switch_time;
	}
}
