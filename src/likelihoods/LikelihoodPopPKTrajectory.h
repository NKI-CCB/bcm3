#pragma once

#include "CVODESolver.h"
#include "ODESolver.h"
#include "Likelihood.h"
#include "RNG.h"
#include "VariableSet.h"
#include "Spinlock.h"

class SBMLModel;

class LikelihoodPopPKTrajectory : public bcm3::Likelihood
{
public:
	LikelihoodPopPKTrajectory(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodPopPKTrajectory();
	
	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize();
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	inline size_t GetNumPatients() const { return patient_ids.size(); }
	inline const VectorReal& GetTimepoints() const { return time; }
	inline const VectorReal& GetObservedConcentrations(size_t patient_id) const { return observed_concentration[patient_id]; }
	inline const VectorReal& GetSimulatedConcentrations(size_t threadix, size_t patient_id) const { return parallel_data[threadix].simulated_concentrations[patient_id]; }
	inline const MatrixReal& GetSimulatedTrajectories(size_t threadix, size_t patient_id) const { return parallel_data[threadix].stored_trajectories[patient_id]; }

private:
	enum PKModelType {
		PKMT_OneCompartment,
		PKMT_TwoCompartment,
		PKMT_TwoCompartmentBiPhasic,
		PKMT_TwoCompartmentLagAndMetabolite,
		PKMT_OneCompartmentTransit,
		PKMT_TwoCompartmentTransit,
		PKMT_Undefined
	};
	struct ParallelData {
		ParallelData();

		Real dose;
		Real dose_cycle2;
		Real dosing_interval;
		unsigned int intermittent;
		std::set<int>* skipped_days;

		Real k_absorption;
		Real k_excretion;
		Real k_elimination;
		Real k_vod;
		Real k_intercompartmental;
		Real k_transit;
		Real k_biphasic_switch_time;
		Real k_absorption2;
		bool biphasic_switch;

		Real current_dose_time;
		MatrixReal simulated_trajectories;
		std::vector<VectorReal> simulated_concentrations;
		std::vector<MatrixReal> stored_trajectories;
	};

	bool CalculateDerivative_OneCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentBiphasic(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_OneCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentLagAndMetabolite(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	Real TreatmentCallback(OdeReal t, void* user);
	Real TreatmentCallbackBiphasic(OdeReal t, void* user);

	// Static variables
	size_t sampling_threads;
	std::shared_ptr<const bcm3::VariableSet> varset;

	VectorReal time;
	std::vector<VectorReal> observed_concentration;
	std::string drug;
	std::vector<Real> dose;
	std::vector<Real> dose_cycle2;
	std::vector<Real> dosing_interval;
	std::vector<unsigned int> intermittent;
	std::vector<std::string> patient_ids;
	std::vector< std::set<int> > skipped_days;
	PKModelType pk_type;
	size_t num_pk_params;
	size_t num_pk_pop_params;

	// Runtime variables
	//std::vector< std::unique_ptr<CVODESolver> > solvers;
	std::vector< std::unique_ptr<ODESolver> > solvers;
	std::vector<ParallelData> parallel_data;

	bcm3::spinlock buffer_spinlock;
	std::vector< std::vector<VectorReal> > prev_parameters;
	std::vector< std::vector<Real> > prev_llh;
	std::vector< size_t > prev_ix;
};
