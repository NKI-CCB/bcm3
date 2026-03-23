#pragma once

#include "Experiment.h"
#include "LinearAlgebraSelector.h"
#include <cvode/cvode.h>

class SBMLModel;

class Cell
{
public:
	Cell(const SBMLModel* model, const Experiment* experiment);
	~Cell();

	bool AllocateSolver(std::string solver_type);
	//bool SetInitialConditionsFromModel(const std::map<size_t, Experiment::SetSpecies>& set_species_map, const std::map<size_t, size_t>& set_init_map, const std::map<size_t, std::vector<int>>& ratio_active_map, const std::map<size_t, std::vector<int>>& ratio_inactive_map, const std::map<size_t, std::vector<size_t>>& ratio_total_active,const std::map<size_t, std::vector<size_t>>& ratio_total_inactive,const VectorReal& transformed_values, Real time);
	bool SetInitialConditionsFromModel();
	bool SetInitialConditionsFromOtherCell(const Cell* other);
	bool Initialize(Real creation_time, const VectorReal& transformed_variables, const VectorReal* sobol_sequence_values, bool is_initial_cell, bool calculate_synchronization_points);

	bool Simulate(Real end_time, Real simulate_past_chromatid_separation_time, VectorReal& output_times, bool &die, bool &divide, Real& achieved_time);

	Real GetInterpolatedSpeciesValue(Real time, size_t i, ESynchronizeCellTrajectory synchronize);
	void RestartInterpolationIteration();
	bool CellAliveAtTime(Real time, ESynchronizeCellTrajectory synchronize) const;
	inline Real GetCreationTime() const { return creation_time; }
	inline bool EnteredMitosis() const { return !std::isnan(nuclear_envelope_breakdown_time); }
	Real GetDuration(EPhaseDuration duration) const;

	size_t GetSolverSteps() const { return solver->GetNumSteps(); }
	OdeReal GetSolverMinStepSize() const { return solver->GetMinStepSize(); }

	static size_t total_num_simulations;
	static size_t cvode_max_steps_reached;
	static size_t cvode_min_timestep_reached;
	static const bool use_generated_code;

private:
	void SetTreatmentConcentration(Real t);
	bool solver_rhs_fn(OdeReal t, const OdeReal* y, OdeReal* ydot, void* user_data);
	bool solver_jac_fn(OdeReal t, const OdeReal* y, const OdeReal* ydot, OdeMatrixReal& jac, void* user_data);
	Real discontinuity_cb(OdeReal t);
	bool integration_step_cb(OdeReal t, const OdeReal* y, void* user_data);

	// Settings
	const SBMLModel* model;
	const Experiment* experiment;
	OdeVectorReal cell_specific_transformed_variables;
	OdeVectorReal cell_specific_non_sampled_transformed_variables;

	bool store_integration_points;
	OdeReal creation_time;				// The time in "experiment time" at which the cell is created. So in "cell time", t=0 will be creation_time
	OdeReal current_simulation_time;

	Real synchronize_offset_time;
	Real simulate_past_chromatid_separation_time;

	// Runtime variables
	size_t DNA_replication_ix;
	size_t DNA_replicated_ix;
	size_t PCNA_gfp_ix;
	size_t nuclear_envelope_ix;
	size_t chromatid_separation_ix;
	size_t cytokinesis_ix;
	size_t apoptosis_ix;

	OdeReal replication_start_time;
	OdeReal replication_finish_time;
	OdeReal PCNA_gfp_increase_time;
	OdeReal nuclear_envelope_breakdown_time;
	OdeReal anaphase_onset_time;

	std::shared_ptr<ODESolver> solver;
	OdeVectorReal initial_y;
	OdeVectorReal constant_species_y;
	OdeVectorReal solver_stored_timepoints;
	OdeMatrixReal solver_output;
	OdeReal simulation_end_time;
	OdeVectorReal simulation_end_y;
	OdeReal previous_integration_step_time;

	bool cell_divided;
	bool cell_died;
};
