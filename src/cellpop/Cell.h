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

	bool SetInitialConditionsFromModel(const std::map<size_t, Experiment::SetSpecies>& set_species_map, Real time);
	bool SetInitialConditionsFromOtherCell(const Cell* other);
	bool Initialize(Real creation_time, const VectorReal& transformed_variables, VectorReal* sobol_sequence_values, bool is_initial_cell, bool calculate_synchronization_point);

	bool Simulate(Real end_time, bool &die, bool &divide, Real& achieved_time);

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
	OdeVectorReal constant_species_y;
	SUNMatrix J;
	bool cvode_initialized;
	Real creation_time;
	OdeReal current_simulation_time;
	size_t cvode_steps;
	Real min_step_size;
	bool completed;

	struct CVodeTimepoint {
		CVodeTimepoint();
		Real cvode_time;
		Real cv_uround;
		Real cv_tn;
		Real cv_h;
		Real cv_hu;
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

	Real replication_start_time;
	Real replication_finish_time;
	Real PCNA_gfp_increase_time;
	Real nuclear_envelope_breakdown_time;
	Real anaphase_onset_time;

	OdeVectorReal cvode_interpolate_y;
	OdeReal interpolation_time;

	derivative_fn derivative;
	jacobian_fn jacobian;
};
