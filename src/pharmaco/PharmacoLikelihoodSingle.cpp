#include "Utils.h"
#include "PharmacoLikelihoodSingle.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include <boost/property_tree/xml_parser.hpp>

PharmacoLikelihoodSingle::PharmacoLikelihoodSingle(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, evaluation_threads(evaluation_threads)
	, additive_sd_ix(std::numeric_limits<size_t>::max())
	, proportional_sd_ix(std::numeric_limits<size_t>::max())
	, use_peripheral_compartment(false)
	, peripheral_forward_rate_ix(std::numeric_limits<size_t>::max())
	, peripheral_backward_rate_ix(std::numeric_limits<size_t>::max())
	, use_transit_compartment(false)
	, num_transit_compartments(0)
	, mean_transit_time_ix(false)
	, dose(std::numeric_limits<Real>::quiet_NaN())
	, dosing_interval(std::numeric_limits<Real>::quiet_NaN())
{
	
}

PharmacoLikelihoodSingle::~PharmacoLikelihoodSingle()
{
}

bool PharmacoLikelihoodSingle::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	bool result = true;

	std::string trial;
	try {
		boost::property_tree::ptree& modelnode = likelihood_node.get_child("pk_model");
		drug = modelnode.get<std::string>("<xmlattr>.drug");
		trial = modelnode.get<std::string>("<xmlattr>.trial");
		patient_id = modelnode.get<std::string>("<xmlattr>.patient", "");
		use_peripheral_compartment = modelnode.get<bool>("<xmlattr>.peripheral_compartment", false);
		num_transit_compartments = modelnode.get<size_t>("<xmlattr>.num_transit_compartments", 0);
	}
	catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	if (num_transit_compartments > 0) {
		use_transit_compartment = true;
	}

	std::string use_patient_str = vm["pk.patient"].as<std::string>();
	if (!use_patient_str.empty()) {
		patient_id = use_patient_str;
	}
	if (patient_id.empty()) {
		LOGERROR("Patient ID has not been specified in either the likelihood or as command-line option.");
		return false;
	}

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

	observation_timepoints.resize(num_timepoints);
	observed_concentrations.resize(num_timepoints);
	result &= data.GetValues(trial, "time", 0, num_timepoints, observation_timepoints);
	result &= data.GetValuesDim2(trial, drug + "_plasma_concentration", patient_ix, 0, num_timepoints, observed_concentrations);
	result &= data.GetValue(trial, drug + "_dose", patient_ix, &dose);
	result &= data.GetValue(trial, drug + "_dosing_interval", patient_ix, &dosing_interval);

	std::vector<Real> times;
	Real last_time = 696;
	Real t = 0;
	while (t < last_time) {
		times.push_back(t);
		t += dosing_interval;
	}
	treatment_timepoints.resize(times.size());
	treatment_doses.resize(times.size());
	for (ptrdiff_t i = 0; i < times.size(); i++) {
		treatment_timepoints(i) = times[i];
		treatment_doses(i) = dose;
	}

	VectorReal new_observed_timepoints = observation_timepoints;
	VectorReal new_observed_concentrations = observed_concentrations;
	ptrdiff_t out_i = 0;
	for (ptrdiff_t i = 0; i < observation_timepoints.size(); i++) {
		if (!std::isnan(observed_concentrations(i))) {
			new_observed_timepoints(out_i) = observation_timepoints(i);
			new_observed_concentrations(out_i) = observed_concentrations(i);
			out_i++;
		}
	}
	observation_timepoints = new_observed_timepoints.segment(0, out_i);
	observed_concentrations = new_observed_concentrations.segment(0, out_i);
	simulated_concentrations.resize(observed_concentrations.size());

	return result;
}

bool PharmacoLikelihoodSingle::PostInitialize()
{
	additive_sd_ix = varset->GetVariableIndex("additive_error_standard_deviation");
	proportional_sd_ix = varset->GetVariableIndex("proportional_error_standard_deviation");
	if (additive_sd_ix == std::numeric_limits<size_t>::max() && proportional_sd_ix == std::numeric_limits<size_t>::max()) {
		LOGERROR("Neither \"additive_error_standard_deviation\" nor \"proportional_error_standard_deviation\" has been specified in the prior; at least one of these variables should be included.");
		return false;
	}

	absorption_ix = varset->GetVariableIndex("absorption");
	excretion_ix = varset->GetVariableIndex("excretion");
	clearance_ix = varset->GetVariableIndex("clearance");
	volume_of_distribution_ix = varset->GetVariableIndex("volume_of_distribution");

	if (absorption_ix == std::numeric_limits<size_t>::max() ||
		excretion_ix == std::numeric_limits<size_t>::max() ||
		clearance_ix == std::numeric_limits<size_t>::max() ||
		volume_of_distribution_ix == std::numeric_limits<size_t>::max()) {
		return false;
	}

	if (use_peripheral_compartment) {
		peripheral_forward_rate_ix = varset->GetVariableIndex("peripheral_forward_rate");
		peripheral_backward_rate_ix = varset->GetVariableIndex("peripheral_backward_rate");
		if (peripheral_forward_rate_ix == std::numeric_limits<size_t>::max() || peripheral_backward_rate_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Peripheral compartmant was specified, but forward or backward rates have not both been specified in prior.");
			return false;
		}
		model.SetUsePeripheralCompartment(true);
	} else {
		model.SetUsePeripheralCompartment(false);
	}

	if (use_transit_compartment) {
		ASSERT(num_transit_compartments > 0);
		mean_transit_time_ix = varset->GetVariableIndex("mean_transit_time");
		if (mean_transit_time_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Transit compartmants were specified, but mean transit time has not been specified in prior.");
			return false;
		}
		model.SetNumTransitCompartments(num_transit_compartments);
	} else {
		model.SetNumTransitCompartments(0);
	}

	return true;
}

bool PharmacoLikelihoodSingle::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	logp = 0.0;

	Real additive_sd = 0.0;
	if (additive_sd_ix != std::numeric_limits<size_t>::max()) {
		additive_sd = varset->TransformVariable(additive_sd_ix, values(additive_sd_ix));
	}
	Real proportional_sd = 0.0;
	if (proportional_sd_ix != std::numeric_limits<size_t>::max()) {
		proportional_sd = varset->TransformVariable(proportional_sd_ix, values(proportional_sd_ix));
	}

	Real absorption = varset->TransformVariable(absorption_ix, values(absorption_ix));
	Real excretion = varset->TransformVariable(excretion_ix, values(excretion_ix));
	Real clearance = varset->TransformVariable(clearance_ix, values(clearance_ix));
	Real volume_of_distribution = varset->TransformVariable(volume_of_distribution_ix, values(volume_of_distribution_ix));

	model.SetAbsorption(absorption);
	model.SetExcretion(excretion);
	model.SetElimination(clearance / volume_of_distribution);

	if (use_peripheral_compartment) {
		Real peripheral_forward = varset->TransformVariable(peripheral_forward_rate_ix, values(peripheral_forward_rate_ix));
		Real peripheral_backward = varset->TransformVariable(peripheral_backward_rate_ix, values(peripheral_backward_rate_ix));
		model.SetPeripheralForwardRate(peripheral_forward);
		model.SetPeripheralBackwardRate(peripheral_backward);
	}
	if (use_transit_compartment) {
		Real mean_transit_time = varset->TransformVariable(mean_transit_time_ix, values(mean_transit_time_ix));
		model.SetTransitRate(num_transit_compartments / mean_transit_time);
	}

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
	concentration_conversion = (1e6 / MW) / volume_of_distribution;

	if (model.Solve(treatment_timepoints, treatment_doses, observation_timepoints, simulated_concentrations, NULL)) {
		for (ptrdiff_t i = 0; i < observation_timepoints.size(); i++) {
			Real x = concentration_conversion * simulated_concentrations(i);
			Real y = observed_concentrations(i);

			if (!std::isnan(y)) {
				logp += bcm3::LogPdfTnu4(x, y, additive_sd + proportional_sd * std::max(x, 0.0));
			}
		}
	} else {
		logp = -std::numeric_limits<Real>::infinity();
	}

	return true;
}

bool PharmacoLikelihoodSingle::GetSimulatedTrajectory(const VectorReal& timepoints, VectorReal& concentrations, MatrixReal& trajectory)
{
	if (!model.Solve(treatment_timepoints, treatment_doses, timepoints, concentrations, &trajectory)) {
		return false;
	}
	for (int i = 0; i < timepoints.size(); i++) {
		concentrations(i) *= concentration_conversion;
	}
	return true;
}

void PharmacoLikelihoodSingle::AddOptionsDescription(boost::program_options::options_description& pod)
{
}
