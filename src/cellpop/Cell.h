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

	bool SetInitialConditionsFromModel(const std::map<size_t, Experiment::SetSpecies>& set_species_map, const std::map<size_t, size_t>& set_init_map, const std::map<size_t, std::vector<int>>& ratio_active_map, const std::map<size_t, std::vector<int>>& ratio_inactive_map, const std::map<size_t, std::vector<size_t>>& ratio_total_active,const std::map<size_t, std::vector<size_t>>& ratio_total_inactive,const VectorReal& transformed_values, Real time);
	bool SetInitialConditionsFromOtherCell(const Cell* other);
	bool Initialize(Real creation_time, const VectorReal& transformed_variables, VectorReal* sobol_sequence_values, bool is_initial_cell, bool calculate_synchronization_point, Real abs_tol, Real rel_tol);

	bool Simulate(Real end_time, Real simulate_past_chromatid_separation_time, bool &die, bool &divide, Real& achieved_time);

	Real GetInterpolatedSpeciesValue(Real time, size_t i, ESynchronizeCellTrajectory synchronize);
	void RestartInterpolationIteration();
	bool CellAliveAtTime(Real time, ESynchronizeCellTrajectory synchronize) const;
	inline bool CellCompleted() const { return completed; }
	inline bool EnteredMitosis() const { return !std::isnan(nuclear_envelope_breakdown_time); }
	inline Real GetCVodeSteps() const { return cvode_steps; }
	inline Real GetCVodeMinStepSize() const { return min_step_size;	}
	Real GetDuration(EPhaseDuration duration) const;

	void SetDerivativeFunctions(derivative_fn fn, jacobian_fn jac);

	static size_t total_num_simulations;
	static size_t cvode_max_steps_reached;
	static size_t cvode_min_timestep_reached;
	static const bool use_generated_code;

private:
	void SetMutations();
	void SetTreatmentConcentration(Real t);
	void RetrieveCVodeInterpolationInfo();
	Real InterpolateEventTime(size_t species_ix, Real threshold, bool above, Real prev_time);
	void CalculateEndY(Real end_time);
	static int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data);
	static int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3);

	const SBMLModel* model;
	const Experiment* experiment;
	OdeVectorReal cell_specific_transformed_variables;
	OdeVectorReal cell_specific_non_sampled_transformed_variables;

	void* cvode_mem;
	SUNLinearSolver LS;
	SUNNonlinearSolver NLS;
	N_Vector cvode_y;
	OdeVectorReal cvode_end_y;
	OdeVectorReal constant_species_y;
	SUNMatrix J;
	bool cvode_initialized;
	OdeReal creation_time;
	OdeReal current_simulation_time;
	size_t cvode_steps;
	OdeReal min_step_size;
	bool completed;

	struct CVodeTimepoint {
		CVodeTimepoint();
		OdeReal cvode_time;
		OdeReal cv_uround;
		OdeReal cv_tn;
		OdeReal cv_h;
		OdeReal cv_hu;
		int cv_q;
	};
	std::vector<CVodeTimepoint> cvode_timepoints;
	OdeMatrixReal cvode_timepoints_zn[6];
	size_t cvode_timepoint_iter;
	Real synchronize_offset_time;

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

	OdeVectorReal cvode_interpolate_y;
	OdeReal interpolation_time;

	derivative_fn derivative;
	jacobian_fn jacobian;
};
