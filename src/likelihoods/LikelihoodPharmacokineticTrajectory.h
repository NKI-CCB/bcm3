#pragma once

#include "ODESolver.h"
#include "Likelihood.h"
#include "RNG.h"
#include "VariableSet.h"

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

	inline const OdeVectorReal& GetTimepoints() const { return time; }
	inline const VectorReal& GetObservedConcentrations() const { return observed_concentration; }
	inline OdeVectorReal GetSimulatedConcentrations(size_t threadix) const { return parallel_data[threadix].simulated_concentrations; }
	inline const OdeMatrixReal& GetSimulatedTrajectories(size_t threadix) const { return parallel_data[threadix].simulated_trajectories; }

	static void AddOptionsDescription(boost::program_options::options_description& pod);

private:
	enum PKModelType {
		PKMT_OneCompartment,
		PKMT_TwoCompartment,
		PKMT_OneCompartmentBiphasicUptake,
		PKMT_TwoCompartmentBiphasicUptake,
		PKMT_OneCompartmentTransit,
		PKMT_TwoCompartmentTransit,
		PKMT_Undefined
	};
	struct ParallelData {
		ParallelData();

		Real dose;
		Real dosing_interval;
		Real dose_after_dose_change;
		Real dose_change_time;
		unsigned int intermittent;
		std::set<int>* skipped_days;

		Real k_absorption;
		Real k_excretion;
		Real k_elimination;
		Real k_vod;
		Real k_periphery_fwd;
		Real k_periphery_bwd;
		Real k_transit;
		Real n_transit;
		Real k_biphasic_switch_time;
		Real k_absorption2;
		bool biphasic_switch;
		Real last_treatment;

		Real current_dose_time;
		OdeMatrixReal simulated_trajectories;
		OdeVectorReal simulated_concentrations;
	};

	bool CalculateDerivative_OneCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartment(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_OneCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_OneCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	bool CalculateDerivative_TwoCompartmentTransit(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);

	bool CalculateJacobian_OneCompartment(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user);
	bool CalculateJacobian_TwoCompartment(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user);
	bool CalculateJacobian_OneCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user);
	bool CalculateJacobian_TwoCompartmentBiphasicUptake(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user);
	bool CalculateJacobian_OneCompartmentTransit(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user);
	bool CalculateJacobian_TwoCompartmentTransit(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeMatrixReal& jac, void* user);

	Real TreatmentCallback(OdeReal t, void* user);
	Real TreatmentCallbackBiphasic(OdeReal t, void* user);
	inline bool CheckGiveTreatment(OdeReal t, ParallelData& pd);
	
	// Static variables
	size_t sampling_threads;
	std::shared_ptr<const bcm3::VariableSet> varset;

	OdeVectorReal time;
	VectorReal observed_concentration;
	std::string drug;
	Real dose;
	Real dosing_interval;
	Real dose_after_dose_change;
	Real dose_change_time;
	bool intermittent;
	std::string patient_id;
	std::set<int> skipped_days;
	PKModelType pk_type;

	Real fixed_vod;
	Real fixed_periphery_fwd;
	Real fixed_periphery_bwd;

	// Runtime variables
	std::vector< std::unique_ptr<ODESolver> > solvers;
	std::vector< ParallelData > parallel_data;
};
