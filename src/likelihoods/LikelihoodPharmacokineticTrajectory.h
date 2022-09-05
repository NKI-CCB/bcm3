#pragma once

#include "CVODESolver.h"
#include "Likelihood.h"
#include "RNG.h"
#include "VariableSet.h"

class SBMLModel;

class LikelihoodPharmacokineticTrajectory : public bcm3::Likelihood
{
public:
	LikelihoodPharmacokineticTrajectory(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodPharmacokineticTrajectory();

	void SetDrugDosing(std::string drug, Real dose, Real interval, bool intermittent);
	void SetPatientID(std::string patient);
	void SetPKModelType(std::string type);
	
	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize();
	virtual bool IsReentrant() { return true; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

private:
	enum PKModelType {
		PKMT_OneCompartment,
		PKMT_TwoCompartment,
		PKMT_TwoCompartmentBiPhasic,
		PKMT_OneCompartmentTransit,
		PKMT_TwoCompartmentTransit,
		PKMT_Undefined
	};
	struct ParallelData {
		Real k_absorption;
		Real k_excretion;
		Real k_elimination;
		Real k_periphery_fwd;
		Real k_periphery_bwd;
		Real k_transit;
		Real k_biphasic_switch_time;
		Real k_absorption2;
		bool biphasic_switch;

		Real current_dose_time;
		MatrixReal simulated_trajectories;
	};

	bool CalculateDerivative_OneCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentBiphasic(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_OneCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	Real TreatmentCallback(OdeReal t, void* user);
	Real TreatmentCallbackBiphasic(OdeReal t, void* user);

	// Static variables
	size_t sampling_threads;
	std::shared_ptr<const bcm3::VariableSet> varset;

	VectorReal time;
	VectorReal observed_concentration;
	std::string drug;
	Real dose;
	Real dosing_interval;
	bool intermittent;
	std::string patient_id;
	PKModelType pk_type;

	// Runtime variables
	std::vector< CVODESolver > solvers;
	std::vector< ParallelData > parallel_data;
};
