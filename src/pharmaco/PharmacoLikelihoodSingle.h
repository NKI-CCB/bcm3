#pragma once

#include "Likelihood.h"
#include "PharmacokineticModel.h"
#include "PharmacoPatient.h"

#include <boost/program_options.hpp>

class PharmacoLikelihoodSingle : public bcm3::Likelihood
{
public:
	PharmacoLikelihoodSingle(size_t sampling_threads, size_t evaluation_threads);
	~PharmacoLikelihoodSingle();
	
	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize();
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	inline const VectorReal& GetTimepoints() const { return patient.observation_timepoints; }
	inline const VectorReal& GetObservedConcentrations() const { return patient.observed_concentrations; }
	inline VectorReal GetSimulatedConcentrations() const { return patient.simulated_concentrations; }
	bool GetSimulatedTrajectory(const VectorReal& timepoints, VectorReal& concentrations, MatrixReal& trajectory);

	static void AddOptionsDescription(boost::program_options::options_description& pod);

private:
	std::shared_ptr<const bcm3::VariableSet> varset;

	std::string drug;
	Real concentration_conversion;

	size_t additive_sd_ix;
	size_t proportional_sd_ix;

	size_t absorption_ix;
	size_t excretion_ix;
	size_t clearance_ix;
	size_t volume_of_distribution_ix;

	PharmacokineticModel model;
	bool use_peripheral_compartment;
	size_t peripheral_forward_rate_ix;
	size_t peripheral_backward_rate_ix;
	bool use_transit_compartment;
	size_t num_transit_compartments;
	size_t mean_transit_time_ix;
	bool biphasic_absorption;
	size_t direct_absorption_ix;
	bool use_metabolite;
	size_t metabolite_conversion_ix;
	size_t metabolite_elimination_ix;

	Patient patient;
};
