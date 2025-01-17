#include "Utils.h"
#include "DrugConstants.h"
#include "PharmacoLikelihoodPopulation.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"

#include <boost/property_tree/xml_parser.hpp>

PharmacoLikelihoodPopulation::ParallelData::ParallelData()
	: concentration_conversion(0.0)
{
}

PharmacoLikelihoodPopulation::PharmacoLikelihoodPopulation(size_t sampling_threads, size_t evaluation_threads)
	: sampling_threads(sampling_threads)
	, use_bioavailability(false)
	, additive_sd_ix(std::numeric_limits<size_t>::max())
	, proportional_sd_ix(std::numeric_limits<size_t>::max())
	, mean_absorption_ix(std::numeric_limits<size_t>::max())
	, mean_excretion_ix(std::numeric_limits<size_t>::max())
	, mean_clearance_ix(std::numeric_limits<size_t>::max())
	, mean_volume_of_distribution_ix(std::numeric_limits<size_t>::max())
	, sigma_absorption_ix(std::numeric_limits<size_t>::max())
	, sigma_excretion_ix(std::numeric_limits<size_t>::max())
	, sigma_clearance_ix(std::numeric_limits<size_t>::max())
	, sigma_volume_of_distribution_ix(std::numeric_limits<size_t>::max())
	, sigma_transit_time_ix(std::numeric_limits<size_t>::max())
	, use_peripheral_compartment(false)
	, peripheral_forward_rate_ix(std::numeric_limits<size_t>::max())
	, peripheral_backward_rate_ix(std::numeric_limits<size_t>::max())
	, use_transit_compartment(false)
	, num_transit_compartments(0)
	, mean_transit_time_ix(std::numeric_limits<size_t>::max())
{
	parallel_data.resize(sampling_threads);
}

PharmacoLikelihoodPopulation::~PharmacoLikelihoodPopulation()
{

}

bool PharmacoLikelihoodPopulation::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	std::string trial;
	try {
		boost::property_tree::ptree& modelnode = likelihood_node.get_child("pk_model");
		drug = modelnode.get<std::string>("<xmlattr>.drug");
		trial = modelnode.get<std::string>("<xmlattr>.trial");
		use_peripheral_compartment = modelnode.get<bool>("<xmlattr>.peripheral_compartment", false);
		num_transit_compartments = modelnode.get<size_t>("<xmlattr>.num_transit_compartments", 0);
		use_bioavailability = modelnode.get<bool>("<xmlattr>.bioavailability", false);
		if (num_transit_compartments > 0) {
			use_transit_compartment = true;
		}
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	bcm3::NetCDFDataFile data;
	if (!data.Open("pkdata.nc", false)) {
		return false;
	}

	bool result = true;
	size_t num_patients;
	std::vector<std::string> patient_ids;
	result &= data.GetDimensionSize(trial, "patients", &num_patients);
	num_patients = 5;
	result &= data.GetValues(trial, "patients", 0, num_patients, patient_ids);

	patients.resize(num_patients);
	for (ptrdiff_t pi = 0; pi < parallel_data.size(); pi++) {
		parallel_data[pi].simulated_concentrations.resize(patients.size());
	}
	for (ptrdiff_t i = 0; i < patients.size(); i++) {
		patients[i].patient_id = patient_ids[i];
		if (!patients[i].Load(data, trial, drug)) {
			return false;
		}
		for (ptrdiff_t pi = 0; pi < parallel_data.size(); pi++) {
			parallel_data[pi].simulated_concentrations[i] = patients[i].simulated_concentrations;
		}
	}

	prev_parameters.resize(parallel_data.size());
	prev_logp.resize(parallel_data.size());
	for (size_t i = 0; i < prev_parameters.size(); i++) {
		prev_parameters[i].resize(num_patients, VectorReal::Zero(10));
		prev_logp[i].resize(num_patients, std::numeric_limits<Real>::quiet_NaN());
	}
	prev_ix.resize(num_patients, 0);

	return true;
}

bool PharmacoLikelihoodPopulation::PostInitialize()
{
	additive_sd_ix = varset->GetVariableIndex("additive_error_standard_deviation", false);
	proportional_sd_ix = varset->GetVariableIndex("proportional_error_standard_deviation", false);
	if (additive_sd_ix == std::numeric_limits<size_t>::max() && proportional_sd_ix == std::numeric_limits<size_t>::max()) {
		LOGERROR("Neither \"additive_error_standard_deviation\" nor \"proportional_error_standard_deviation\" has been specified in the prior; at least one of these variables should be included.");
		return false;
	}

	for (ptrdiff_t i = 0; i < parallel_data.size(); i++) {
		parallel_data[i].model.SetUsePeripheralCompartment(use_peripheral_compartment);
		if (use_transit_compartment) {
			parallel_data[i].model.SetNumTransitCompartments(num_transit_compartments);
		} else {
			parallel_data[i].model.SetNumTransitCompartments(0);
		}
		parallel_data[i].cache_lookup_params.setZero(10);
	}

	mean_absorption_ix = varset->GetVariableIndex("mean_absorption");
	mean_clearance_ix = varset->GetVariableIndex("mean_clearance");
	mean_volume_of_distribution_ix = varset->GetVariableIndex("mean_volume_of_distribution");

	if (mean_absorption_ix == std::numeric_limits<size_t>::max() ||
		mean_clearance_ix == std::numeric_limits<size_t>::max() ||
		mean_volume_of_distribution_ix == std::numeric_limits<size_t>::max()) {
		return false;
	}

	mean_excretion_ix = varset->GetVariableIndex("mean_excretion", false);
	sigma_absorption_ix = varset->GetVariableIndex("sigma_absorption", false);
	sigma_excretion_ix = varset->GetVariableIndex("sigma_excretion", false);
	sigma_clearance_ix = varset->GetVariableIndex("sigma_clearance", false);
	sigma_volume_of_distribution_ix = varset->GetVariableIndex("sigma_volume_of_distribution", false);
	sigma_transit_time_ix = varset->GetVariableIndex("sigma_transit_time", false);

	if (sigma_absorption_ix != std::numeric_limits<size_t>::max()) {
		if (!InitializePatientMarginals("absorption", patient_absorption_ix)) {
			return false;
		}
	}
	if (sigma_excretion_ix != std::numeric_limits<size_t>::max()) {
		if (!InitializePatientMarginals("excretion", patient_excretion_ix)) {
			return false;
		}
	}
	if (sigma_clearance_ix != std::numeric_limits<size_t>::max()) {
		if (!InitializePatientMarginals("clearance", patient_clearance_ix)) {
			return false;
		}
	}
	if (sigma_volume_of_distribution_ix != std::numeric_limits<size_t>::max()) {
		if (!InitializePatientMarginals("volume_of_distribution", patient_volume_of_distribution_ix)) {
			return false;
		}
	}
	if (sigma_transit_time_ix != std::numeric_limits<size_t>::max()) {
		if (!InitializePatientMarginals("transit_time", patient_transit_time_ix)) {
			return false;
		}
	}

	if (use_peripheral_compartment) {
		peripheral_forward_rate_ix = varset->GetVariableIndex("peripheral_forward_rate");
		peripheral_backward_rate_ix = varset->GetVariableIndex("peripheral_backward_rate");
		if (peripheral_forward_rate_ix == std::numeric_limits<size_t>::max() || peripheral_backward_rate_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Peripheral compartmant was specified, but forward or backward rates have not both been specified in prior.");
			return false;
		}
	}

	if (use_transit_compartment) {
		ASSERT(num_transit_compartments > 0);
		mean_transit_time_ix = varset->GetVariableIndex("mean_transit_time");
		if (mean_transit_time_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Transit compartmants were specified, but mean transit time has not been specified in prior.");
			return false;
		}
	}

	if (use_bioavailability) {
		if (!InitializePatientMarginals("bioavailability", patient_bioavailability_ix)) {
			return false;
		}
	}

	for (size_t i = 0; i < parallel_data.size(); i++) {
		parallel_data[i].model.SetUsePeripheralCompartment(use_peripheral_compartment);
		parallel_data[i].model.SetNumTransitCompartments(num_transit_compartments);
	}

	return true;
}

bool PharmacoLikelihoodPopulation::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	ParallelData& pd = parallel_data[threadix];

	logp = 0.0;

	Real additive_sd = 1e-12;
	if (additive_sd_ix != std::numeric_limits<size_t>::max()) {
		additive_sd = varset->TransformVariable(additive_sd_ix, values(additive_sd_ix));
	}
	Real proportional_sd = 0.0;
	if (proportional_sd_ix != std::numeric_limits<size_t>::max()) {
		proportional_sd = varset->TransformVariable(proportional_sd_ix, values(proportional_sd_ix));
	}

	for (ptrdiff_t i = 0; i < patients.size(); i++) {
		Real this_logp = 0.0;
		const Patient& patient = patients[i];

		SetupSimulation(threadix, values, i);

		pd.cache_lookup_params(8) = additive_sd;
		pd.cache_lookup_params(9) = proportional_sd;
		if (!LookupCache(pd.cache_lookup_params, i, this_logp)) {
			bool result = pd.model.Solve(patient.treatment_timepoints, patient.treatment_doses, patient.observation_timepoints, pd.simulated_concentrations[i], nullptr);
			if (result) {
				for (ptrdiff_t ti = 0; ti < patient.observation_timepoints.size(); ti++) {
					Real x = pd.concentration_conversion * pd.simulated_concentrations[i](ti);
					if (std::isnan(x)) {
						this_logp = -std::numeric_limits<Real>::infinity();
						break;
					}

					Real y = patient.observed_concentrations(ti);
					if (!std::isnan(y)) {
						this_logp += bcm3::LogPdfTnu4(x, y, additive_sd + proportional_sd * std::max(x, 0.0));
					}
				}
			} else {
				this_logp = -std::numeric_limits<Real>::infinity();
				break;
			}

			SetCache(pd.cache_lookup_params, i, this_logp);
		}
		logp += this_logp;
	}

	return true;
}

bool PharmacoLikelihoodPopulation::GetSimulatedTrajectory(size_t threadix, const VectorReal& values, size_t patient_ix, const VectorReal& timepoints, VectorReal& concentrations, MatrixReal& trajectory)
{
	if (patient_ix > patients.size()) {
		return false;
	}

	ParallelData& pd = parallel_data[threadix];
	const Patient& patient = patients[patient_ix];

	SetupSimulation(threadix, values, patient_ix);
	bool result = pd.model.Solve(patient.treatment_timepoints, patient.treatment_doses, timepoints, concentrations, &trajectory);
	if (result) {
		concentrations *= pd.concentration_conversion;
		return true;
	} else {
		return false;
	}
}

void PharmacoLikelihoodPopulation::SetupSimulation(size_t threadix, const VectorReal& values, size_t patient_ix)
{
	ParallelData& pd = parallel_data[threadix];
	const Patient& patient = patients[patient_ix];

	Real absorption;
	if (sigma_absorption_ix == std::numeric_limits<size_t>::max()) {
		absorption = bcm3::fastpow10(values(mean_absorption_ix));
	} else {
		absorption = bcm3::fastpow10(bcm3::QuantileNormal(values(patient_absorption_ix[patient_ix]), values(mean_absorption_ix), values(sigma_absorption_ix)));
	}

	Real excretion;
	if (mean_excretion_ix == std::numeric_limits<size_t>::max()) {
		excretion = 0.0;
	} else {
		if (sigma_excretion_ix == std::numeric_limits<size_t>::max()) {
			excretion = bcm3::fastpow10(values(mean_excretion_ix));
		} else {
			excretion = bcm3::fastpow10(bcm3::QuantileNormal(values(patient_excretion_ix[patient_ix]), values(mean_excretion_ix), values(sigma_excretion_ix)));
		}
	}

	Real clearance;
	if (sigma_clearance_ix == std::numeric_limits<size_t>::max()) {
		clearance = bcm3::fastpow10(values(mean_clearance_ix));
	} else {
		clearance = bcm3::fastpow10(bcm3::QuantileNormal(values(patient_clearance_ix[patient_ix]), values(mean_clearance_ix), values(sigma_clearance_ix)));
	}

	Real volume_of_distribution;
	if (sigma_volume_of_distribution_ix == std::numeric_limits<size_t>::max()) {
		volume_of_distribution = bcm3::fastpow10(values(mean_volume_of_distribution_ix));
	} else {
		volume_of_distribution = bcm3::fastpow10(bcm3::QuantileNormal(values(patient_volume_of_distribution_ix[patient_ix]), values(mean_volume_of_distribution_ix), values(sigma_volume_of_distribution_ix)));
	}

	pd.cache_lookup_params(0) = absorption;
	pd.cache_lookup_params(1) = excretion;
	pd.cache_lookup_params(2) = clearance;
	pd.cache_lookup_params(3) = volume_of_distribution;

	if (use_peripheral_compartment) {
		Real peripheral_forward = varset->TransformVariable(peripheral_forward_rate_ix, values(peripheral_forward_rate_ix));
		Real peripheral_backward = varset->TransformVariable(peripheral_backward_rate_ix, values(peripheral_backward_rate_ix));
		pd.model.SetPeripheralForwardRate(peripheral_forward);
		pd.model.SetPeripheralBackwardRate(peripheral_backward);
		pd.cache_lookup_params(4) = peripheral_forward;
		pd.cache_lookup_params(5) = peripheral_backward;
	}
	if (use_transit_compartment) {
		Real transit_time;
		if (sigma_transit_time_ix == std::numeric_limits<size_t>::max()) {
			transit_time = varset->TransformVariable(mean_transit_time_ix, values(mean_transit_time_ix));
		} else {
			transit_time = bcm3::fastpow10(bcm3::QuantileNormal(values(patient_transit_time_ix[patient_ix]), values(mean_transit_time_ix), values(sigma_transit_time_ix)));
		}
		pd.model.SetTransitRate(num_transit_compartments / transit_time);
		pd.cache_lookup_params(6) = transit_time;
	}
	if (use_bioavailability) {
		Real bioavailability = values(patient_bioavailability_ix[patient_ix]);
		pd.model.SetBioavailability(bioavailability);
		pd.cache_lookup_params(7) = bioavailability;
	}
	pd.model.SetAbsorption(absorption);
	pd.model.SetExcretion(excretion);
	pd.model.SetElimination(clearance / volume_of_distribution);
	pd.concentration_conversion = (1e6 / GetDrugMolecularWeight(drug)) / volume_of_distribution;
}

bool PharmacoLikelihoodPopulation::InitializePatientMarginals(std::string name, std::vector<size_t>& ixs)
{
	ixs.resize(patients.size());
	for (ptrdiff_t i = 0; i < patients.size(); i++) {
		std::string varname = std::string("p") + std::to_string(i + 1) + "_" + name;
		ixs[i] = varset->GetVariableIndex(varname);
		if (ixs[i] == std::numeric_limits<size_t>::max()) {
			LOGERROR("Standard deviation found for \"%s\", but could not find prior variable for \"%s\"", name.c_str(), varname.c_str());
			return false;
		}
	}
	return true;
}

bool PharmacoLikelihoodPopulation::LookupCache(const VectorReal& params, size_t patient_ix, Real& logp)
{
	buffer_spinlock.lock();
	ptrdiff_t which_exactly_equal = -1;
	for (ptrdiff_t hi = 0; hi < prev_parameters.size(); hi++) {
		bool exactly_equal = true;
		for (int i = 0; i < params.size(); i++) {
			if (params(i) != prev_parameters[hi][patient_ix](i)) {
				exactly_equal = false;
				break;
			}
		}
		if (exactly_equal) {
			which_exactly_equal = hi;
			break;
		}
	}
	if (which_exactly_equal != -1) {
		logp = prev_logp[which_exactly_equal][patient_ix];
		buffer_spinlock.unlock();
		return true;
	} else {
		buffer_spinlock.unlock();
		return false;
	}
}

void PharmacoLikelihoodPopulation::SetCache(const VectorReal& params, size_t patient_ix, Real logp)
{
	buffer_spinlock.lock();
	prev_parameters[prev_ix[patient_ix]][patient_ix] = params;
	prev_logp[prev_ix[patient_ix]][patient_ix] = logp;
	prev_ix[patient_ix]++;
	if (prev_ix[patient_ix] >= prev_parameters.size()) {
		prev_ix[patient_ix] = 0;
	}
	buffer_spinlock.unlock();
}
