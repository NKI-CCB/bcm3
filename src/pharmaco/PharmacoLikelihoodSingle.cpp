#include "Utils.h"
#include "DrugConstants.h"
#include "PharmacoLikelihoodSingle.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include <boost/property_tree/xml_parser.hpp>

PharmacoLikelihoodSingle::PharmacoLikelihoodSingle(size_t sampling_threads, size_t evaluation_threads)
	: concentration_conversion(0.0)
	, additive_sd_ix(std::numeric_limits<size_t>::max())
	, proportional_sd_ix(std::numeric_limits<size_t>::max())
	, absorption_ix(std::numeric_limits<size_t>::max())
	, excretion_ix(std::numeric_limits<size_t>::max())
	, clearance_ix(std::numeric_limits<size_t>::max())
	, volume_of_distribution_ix(std::numeric_limits<size_t>::max())
	, use_peripheral_compartment(false)
	, peripheral_forward_rate_ix(std::numeric_limits<size_t>::max())
	, peripheral_backward_rate_ix(std::numeric_limits<size_t>::max())
	, use_transit_compartment(false)
	, num_transit_compartments(0)
	, mean_transit_time_ix(std::numeric_limits<size_t>::max())
{
	
}

PharmacoLikelihoodSingle::~PharmacoLikelihoodSingle()
{
}

bool PharmacoLikelihoodSingle::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	std::string trial;
	try {
		boost::property_tree::ptree& modelnode = likelihood_node.get_child("pk_model");
		drug = modelnode.get<std::string>("<xmlattr>.drug");
		trial = modelnode.get<std::string>("<xmlattr>.trial");
		patient.patient_id = modelnode.get<std::string>("<xmlattr>.patient", "");
		use_peripheral_compartment = modelnode.get<bool>("<xmlattr>.peripheral_compartment", false);
		num_transit_compartments = modelnode.get<size_t>("<xmlattr>.num_transit_compartments", 0);
		if (num_transit_compartments > 0) {
			use_transit_compartment = true;
		}
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	std::string use_patient_str = vm["pharmacosingle.patient"].as<std::string>();
	if (!use_patient_str.empty()) {
		patient.patient_id = use_patient_str;
	}
	if (patient.patient_id.empty()) {
		LOGERROR("Patient ID has not been specified in either the likelihood or as command-line option.");
		return false;
	}

	bcm3::NetCDFDataFile data;
	if (!data.Open("pkdata.nc", false)) {
		return false;
	}
	if (!patient.Load(data, trial, drug)) {
		return false;
	}

	return true;
}

bool PharmacoLikelihoodSingle::PostInitialize()
{
	additive_sd_ix = varset->GetVariableIndex("additive_error_standard_deviation", false);
	proportional_sd_ix = varset->GetVariableIndex("proportional_error_standard_deviation", false);
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

	concentration_conversion = (1e6 / GetDrugMolecularWeight(drug)) / volume_of_distribution;

	if (model.Solve(patient.treatment_timepoints, patient.treatment_doses, patient.observation_timepoints, patient.simulated_concentrations, NULL)) {
		for (ptrdiff_t i = 0; i < patient.observation_timepoints.size(); i++) {
			Real x = concentration_conversion * patient.simulated_concentrations(i);
			Real y = patient.observed_concentrations(i);

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
	if (!model.Solve(patient.treatment_timepoints, patient.treatment_doses, timepoints, concentrations, &trajectory)) {
		return false;
	}
	for (int i = 0; i < timepoints.size(); i++) {
		concentrations(i) *= concentration_conversion;
	}
	return true;
}

void PharmacoLikelihoodSingle::AddOptionsDescription(boost::program_options::options_description& pod)
{
	pod.add_options()
		("pharmacosingle.patient", boost::program_options::value<std::string>()->default_value(""), "Fit the data from this specific patient.")
	;
}
