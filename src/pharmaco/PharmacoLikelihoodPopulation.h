#pragma once

#include "Likelihood.h"
#include "PharmacokineticModel.h"
#include "PharmacoPatient.h"
#include "Spinlock.h"

#include <boost/program_options.hpp>

class PharmacoLikelihoodPopulation : public bcm3::Likelihood
{
public:
	PharmacoLikelihoodPopulation(size_t sampling_threads, size_t evaluation_threads);
	~PharmacoLikelihoodPopulation();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize();
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	inline const VectorReal& GetTimepoints(size_t patient_ix) const { return patients[patient_ix].observation_timepoints; }
	inline const VectorReal& GetObservedConcentrations(size_t patient_ix) const { return patients[patient_ix].observed_concentrations; }
	inline size_t GetNumPatients() const { return patients.size(); }
	bool GetSimulatedTrajectory(size_t threadix, const VectorReal& values, size_t patient_ix, const VectorReal& timepoints, VectorReal& concentrations, MatrixReal& trajectory); // This is not reentrant and assumes that EvaluateLogProbability has just been called, from threadix 0

private:
	void SetupSimulation(size_t threadix, const VectorReal& values, size_t patient_ix);
	bool InitializePatientMarginals(std::string name, std::vector<size_t>& ixs);
	bool LookupCache(const VectorReal& params, size_t patient_ix, Real& logp);
	void SetCache(const VectorReal& params, size_t patient_ix, Real logp);

	struct ParallelData
	{
		ParallelData();
		PharmacokineticModel model;
		Real concentration_conversion;
		std::vector<VectorReal> simulated_concentrations;
		VectorReal cache_lookup_params;
	};

	size_t sampling_threads;
	std::shared_ptr<const bcm3::VariableSet> varset;

	std::string drug;
	bool use_bioavailability;

	size_t additive_sd_ix;
	size_t proportional_sd_ix;

	size_t mean_absorption_ix;
	size_t mean_excretion_ix;
	size_t mean_clearance_ix;
	size_t mean_volume_of_distribution_ix;
	size_t sigma_absorption_ix;
	size_t sigma_excretion_ix;
	size_t sigma_clearance_ix;
	size_t sigma_volume_of_distribution_ix;
	size_t sigma_transit_time_ix;

	std::vector<size_t> patient_absorption_ix;
	std::vector<size_t> patient_excretion_ix;
	std::vector<size_t> patient_clearance_ix;
	std::vector<size_t> patient_volume_of_distribution_ix;
	std::vector<size_t> patient_bioavailability_ix;
	std::vector<size_t> patient_transit_time_ix;

	bool use_peripheral_compartment;
	size_t peripheral_forward_rate_ix;
	size_t peripheral_backward_rate_ix;
	bool use_transit_compartment;
	size_t num_transit_compartments;
	size_t mean_transit_time_ix;

	std::vector<ParallelData> parallel_data;
	std::vector<Patient> patients;

	bcm3::spinlock buffer_spinlock;
	std::vector< std::vector<VectorReal> > prev_parameters;
	std::vector< std::vector<Real> > prev_logp;
	std::vector< size_t > prev_ix;
};
