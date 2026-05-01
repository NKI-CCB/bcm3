#include "Utils.h"
#include "Cell.h"
#include "ProbabilityDistributions.h"
#include "SBMLAssignmentRule.h"
#include "SBMLModel.h"
#include "SBMLSpecies.h"
#include "VariabilityDescription.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <fstream>
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"
#include "ODESolverDP5.h"

size_t Cell::total_num_simulations = 0;
size_t Cell::cvode_max_steps_reached = 0;
size_t Cell::cvode_min_timestep_reached = 0;
const bool Cell::use_generated_code = 1;
static const int max_cvode_steps = 10000;

Cell::Cell(const SBMLModel* model, const Experiment* experiment)
	: model(model)
	, experiment(experiment)
	, store_integration_points(false)
	, creation_time(std::numeric_limits<Real>::quiet_NaN())
	, current_simulation_time(std::numeric_limits<Real>::quiet_NaN())
	, synchronize_offset_time(0.0)
	, DNA_replication_ix(std::numeric_limits<size_t>::max())
	, DNA_replicated_ix(std::numeric_limits<size_t>::max())
	, PCNA_gfp_ix(std::numeric_limits<size_t>::max())
	, nuclear_envelope_ix(std::numeric_limits<size_t>::max())
	, chromatid_separation_ix(std::numeric_limits<size_t>::max())
	, cytokinesis_ix(std::numeric_limits<size_t>::max())
	, apoptosis_ix(std::numeric_limits<size_t>::max())
	, replication_start_time(std::numeric_limits<Real>::quiet_NaN())
	, replication_finish_time(std::numeric_limits<Real>::quiet_NaN())
	, PCNA_gfp_increase_time(std::numeric_limits<Real>::quiet_NaN())
	, nuclear_envelope_breakdown_time(std::numeric_limits<Real>::quiet_NaN())
	, anaphase_onset_time(std::numeric_limits<Real>::quiet_NaN())
{
	initial_y.setZero(model->GetNumODEIntegratedSpecies());
	constant_species_y.setZero(model->GetNumConstantSpecies());
	cell_specific_transformed_variables.setConstant(experiment->GetVarset()->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());
	cell_specific_non_sampled_transformed_variables.setConstant(experiment->non_sampled_parameters.size(), std::numeric_limits<Real>::quiet_NaN());

	DNA_replication_ix = model->GetSimulatedSpeciesByName("replicating_DNA", false);
	DNA_replicated_ix = model->GetSimulatedSpeciesByName("replicated_DNA", false);
	PCNA_gfp_ix = model->GetSimulatedSpeciesByName("PCNA_gfp", false);
	nuclear_envelope_ix = model->GetSimulatedSpeciesByName("nuclear_envelope", false);
	chromatid_separation_ix = model->GetSimulatedSpeciesByName("chromatid_separation", false);
	cytokinesis_ix = model->GetSimulatedSpeciesByName("cytokinesis", false);
	apoptosis_ix = model->GetSimulatedSpeciesByName("apoptosis", false);
}

Cell::~Cell()
{
}

bool Cell::AllocateSolver(std::string solver_type)
{
	if (solver_type == "DP5") {
		solver = std::make_shared<ODESolverDP5>();
	} else if (solver_type == "CVODE") {
		solver = std::make_shared<ODESolverCVODE>();
	} else {
		LOGERROR("Unknown solver type \"%s\"; accepted options are \"DP5\" or \"CVODE\"", solver_type.c_str());
		return false;
	}

	solver->SetDerivativeFunction(boost::bind(&Cell::solver_rhs_fn, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4));
	solver->SetIntegrationStepCallback(boost::bind(&Cell::integration_step_cb, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3));
	solver->Initialize(model->GetNumODEIntegratedSpecies(), (void*)this, experiment->solver_max_steps);
	solver->SetTolerance(experiment->solver_rel_tol, experiment->solver_abs_tol);
	solver->SetSolverParameter("min_dt", 0, experiment->solver_min_timestep);
	solver->SetSolverParameter("max_dt", 0, experiment->solver_max_timestep);
	solver->SetSolverParameter("max_steps", experiment->solver_max_steps, std::numeric_limits<OdeReal>::quiet_NaN());
	return true;
}

//bool Cell::SetInitialConditionsFromModel(const std::map<size_t, Experiment::SetSpecies>& set_species_map, const std::map<size_t, size_t>& set_init_map, const std::map<size_t, std::vector<int>>& ratio_active_map,const std::map<size_t, std::vector<int>>& ratio_inactive_map, const std::map<size_t, std::vector<size_t>>& ratio_total_active, const std::map<size_t, std::vector<size_t>>& ratio_total_inactive,const VectorReal& transformed_values, Real time)
bool Cell::SetInitialConditionsFromModel()
{
	for (size_t i = 0; i < model->GetNumODEIntegratedSpecies(); i++) {
		initial_y(i) = model->GetODEIntegratedSpecies(i)->GetInitialValue();
	}
	for (size_t i = 0; i < model->GetNumConstantSpecies(); i++) {
		constant_species_y(i) = model->GetConstantSpecies(i)->GetInitialValue();
	}

#if 0
	for (std::map<size_t, Experiment::SetSpecies>::const_iterator ssmi = set_species_map.begin(); ssmi != set_species_map.end(); ++ssmi) {
		if (time >= ssmi->second.begin_time && time < ssmi->second.end_time) {
			initial_y(ssmi->first) = ssmi->second.value;
		}
	}

	for(auto const& sic : set_init_map){
		initial_y(sic.first) = transformed_values[sic.second];
	}

	for(auto const& ram : ratio_active_map){
		initial_y(ram.first) = transformed_values[ram.second[0]] * transformed_values[ram.second[1]];
	}

	for(auto const& rim : ratio_inactive_map){
		initial_y(rim.first) = (1 - transformed_values[rim.second[0]]) * transformed_values[rim.second[1]];
	}

	for(auto const& rtm : ratio_total_active){
		initial_y(rtm.first) = transformed_values[rtm.second[0]] * (model -> GetODEIntegratedSpecies(rtm.second[1]) -> GetInitialValue() + model ->GetODEIntegratedSpecies(rtm.second[2]) -> GetInitialValue());
	}

	for(auto const& rtm : ratio_total_inactive){
		initial_y(rtm.first) = (1 - transformed_values[rtm.second[0]]) * (model ->GetODEIntegratedSpecies(rtm.second[1]) -> GetInitialValue() + model ->GetODEIntegratedSpecies(rtm.second[2]) -> GetInitialValue());
	}
#endif

	return true;
}

bool Cell::SetInitialConditionsFromOtherCell(const Cell* other)
{
	for (size_t i = 0; i < model->GetNumODEIntegratedSpecies(); i++) {
		initial_y(i) = other->simulation_end_y(i);
	}
	constant_species_y = other->constant_species_y;

#if 0
	initial_y(model->GetODEIntegratedSpecies("DNA", true)) *= 0.5;
	initial_y(model->GetODEIntegratedSpecies("replicating_DNA", true)) *= 0.5;
	initial_y(model->GetODEIntegratedSpecies("replicated_DNA", true)) *= 0.5;
	initial_y(model->GetODEIntegratedSpecies("licensed_DNA", true)) *= 0.5;
	initial_y(model->GetODEIntegratedSpecies("mitogenic_signal", true)) = 0.0;
	initial_y(model->GetODEIntegratedSpecies("unassembled_spindle", true)) = 0.0;
	initial_y(model->GetODEIntegratedSpecies("assembled_spindle", true)) = 0.0;
	initial_y(model->GetODEIntegratedSpecies("chromatid_separation", true)) = 0.0;
	initial_y(model->GetODEIntegratedSpecies("nuclear_envelope", true)) = 1.0;
	initial_y(model->GetODEIntegratedSpecies("cytokinesis", true)) = 0.0;
	initial_y(model->GetODEIntegratedSpecies("G2_delay", true)) = 0.0;
#endif

	return true;
}

bool Cell::Initialize(Real creation_time, const VectorReal& transformed_variables, const VectorReal* sobol_sequence_values, bool is_initial_cell, bool calculate_synchronization_points)
{
	store_integration_points = calculate_synchronization_points;

	// Make a copy of the sampled and non-sampled parameters so that we can apply variability to them
	for (size_t i = 0; i < cell_specific_transformed_variables.size(); i++) {
		cell_specific_transformed_variables(i) = transformed_variables(i);
	}
	for (size_t i = 0; i < cell_specific_non_sampled_transformed_variables.size(); i++) {
		cell_specific_non_sampled_transformed_variables(i) = experiment->non_sampled_parameters(i);
	}

	// Apply variabilities
	int sobol_sequence_ix = 0;
	for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
		VectorReal pseudorandom_vector = (*it)->GetPseudorandomVector(*sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters);

		for (size_t i = 0; i < cell_specific_transformed_variables.size(); i++) {
			(*it)->ApplyVariabilityParameter(experiment->varset->GetVariableName(i), cell_specific_transformed_variables(i), pseudorandom_vector, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
		}
		for (size_t i = 0; i < cell_specific_non_sampled_transformed_variables.size(); i++) {
			(*it)->ApplyVariabilityParameter(experiment->non_sampled_parameter_names[i], cell_specific_non_sampled_transformed_variables(i), pseudorandom_vector, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
		}
		for (size_t i = 0; i < model->GetNumODEIntegratedSpecies(); i++) {
			(*it)->ApplyVariabilityInitialCondition(model->GetODEIntegratedSpeciesName(i), initial_y(i), pseudorandom_vector, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
		}
	}

	this->creation_time = creation_time;

	replication_start_time = std::numeric_limits<Real>::quiet_NaN();
	replication_finish_time = std::numeric_limits<Real>::quiet_NaN();
	PCNA_gfp_increase_time = std::numeric_limits<Real>::quiet_NaN();
	nuclear_envelope_breakdown_time = std::numeric_limits<Real>::quiet_NaN();
	anaphase_onset_time = std::numeric_limits<Real>::quiet_NaN();

	total_num_simulations++;

	solver->Restart();

	return true;
}

bool Cell::Simulate(Real end_time, Real simulate_past_chromatid_separation_time, VectorReal& output_times, bool& die, bool& divide, Real& achieved_time)
{
	bool result = true;

	cell_died = false;
	cell_divided = false;
	current_simulation_time = 0;
	previous_integration_step_time = 0;
	simulation_end_time = (OdeReal)end_time - creation_time;
	this->simulate_past_chromatid_separation_time = simulate_past_chromatid_separation_time;

	if (output_times.size() > 0) {
		solver_stored_timepoints = OdeVectorReal(output_times.size());
		for (int i = 0; i < output_times.size(); i++) {
			solver_stored_timepoints(i) = output_times(i) - creation_time;
		}
	}
	
	// Find any discontinuities in the treatment trajectories, and inform the solver of the first one
	Real first_discontinuity = std::numeric_limits<Real>::quiet_NaN();
	for (int i = 0; i < experiment->treatment_trajectories.size(); i++) {
		Real discontinuity = experiment->treatment_trajectories[i]->FirstDiscontinuity(creation_time);
		if (!std::isnan(discontinuity)) {
			while (discontinuity < 0.0) {
				discontinuity = experiment->treatment_trajectories[i]->NextDiscontinuity(discontinuity, creation_time);
			}
		}

		if (!(first_discontinuity < discontinuity)) { // Slightly odd comparison to handle nan's
			first_discontinuity = discontinuity;
		}
	}
	if (!std::isnan(first_discontinuity) && first_discontinuity > 0.0) {
		ODESolver::TDiscontinuityCallback cb = boost::bind(&Cell::discontinuity_cb, this, boost::placeholders::_1);
		solver->SetDiscontinuity(first_discontinuity, cb, nullptr);
	}

	// Integrate the ODE system
	if (store_integration_points) {
		result = solver->SolveStoreIntegrationPoints(initial_y, output_times.tail<1>()[0]);
	} else {
		result = solver->SolveReturnSolution(initial_y, &solver_stored_timepoints, &solver_output);
		
	}

	divide = cell_divided;
	die = cell_died;

	Real achieved_cell_time;
	if (cell_divided || cell_died) {
		// Simulation ended prematurely due to division or death
		achieved_cell_time = simulation_end_time;
	} else {
		if (store_integration_points) {
			simulation_end_y = solver->get_current_y();
			achieved_cell_time = previous_integration_step_time;
		} else {
			simulation_end_y = solver_output.col(solver_stored_timepoints.size() - 1);
			achieved_cell_time = solver_stored_timepoints.tail<1>()[0];
		}
	}

	if (result) {
#if 0
		for (int ti = 0; ti < solver_stored_timepoints.size(); ti++) {
			static_assert(!solver_output.IsRowMajor);
			OdeReal* y = solver_output.col(ti).data();
			for (size_t i = 0; i < model->GetNumAssignmentRules(); i++) {
				OdeReal value;
				model->GetAssignmentRule(i).Calculate(y, nullptr, cell_specific_transformed_variables.data(), cell_specific_non_sampled_transformed_variables.data(), &value);

				constant_solver_output
			}
		}
#endif
	}
	achieved_time = achieved_cell_time + creation_time;

	return result;
}

void Cell::RestartInterpolationIteration()
{
	solver->RestartInterpolationIteration();
}

Real Cell::GetInterpolatedSpeciesValue(Real time, size_t species_ix, ESynchronizeCellTrajectory synchronize)
{
	if (species_ix >= model->GetNumSimulatedSpecies()) {
		LOGERROR("Out of bounds species index");
		ASSERT(false);
		return std::numeric_limits<Real>::quiet_NaN();
	}
	if (synchronize != ESynchronizeCellTrajectory::None) {
		// Can only synchronize if integration points were stored
		ASSERT(store_integration_points);
	}

	Real cell_time;
	if (synchronize == ESynchronizeCellTrajectory::DNAReplicationStart) {
		if (std::isnan(replication_start_time)) {
			cell_time = time + simulation_end_time;
		} else {
			cell_time = time + replication_start_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::PCNA_gfp_increase) {
		if (std::isnan(PCNA_gfp_increase_time)) {
			cell_time = time + simulation_end_time;
		} else {
			cell_time = time + PCNA_gfp_increase_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::NuclearEnvelopeBreakdown) {
		if (std::isnan(nuclear_envelope_breakdown_time)) {
			cell_time = time + simulation_end_time;
		} else {
			cell_time = time + nuclear_envelope_breakdown_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::AnaphaseOnset) {
		if (std::isnan(anaphase_onset_time)) {
			cell_time = time + simulation_end_time;
		} else {
			cell_time = time + anaphase_onset_time;
		}
	} else {
		cell_time = time - creation_time;
	}
	if (cell_time < 0.0) {
		return std::numeric_limits<Real>::quiet_NaN();
	}

	if (species_ix < model->GetNumODEIntegratedSpecies()) {
		if (store_integration_points) {
			// Integration time points were stored, use them to interpolate
			return solver->GetInterpolatedY(cell_time, species_ix);
		} else {
			// Exact timepoints should be available in the solver output
			for (int i = 0; i < solver_stored_timepoints.size(); i++) {
				if (solver_stored_timepoints(i) == cell_time) {
					return solver_output(species_ix, i);
				}
			}
			return std::numeric_limits<Real>::quiet_NaN();
		}
	} else {
		SetTreatmentConcentration(cell_time);

		if (store_integration_points) {
			// Integration time points were stored, use them to interpolate
			return solver->GetInterpolatedY(cell_time, species_ix);
		} else {
			// Exact timepoints should be available in the solver output
			for (int i = 0; i < solver_stored_timepoints.size(); i++) {
				if (solver_stored_timepoints(i) == cell_time) {
					size_t constant_species_ix = species_ix - model->GetNumODEIntegratedSpecies();
					size_t simulated_species_ix = model->GetSimulatedSpeciesFromConstantSpecies(constant_species_ix);
					const SBMLAssignmentRule* assignment_rule = model->GetAssignmentRuleForSimulatedSpecies(simulated_species_ix);
					if (assignment_rule) {
						OdeReal value;
						assignment_rule->Calculate(solver_output.col(i).data(), constant_species_y.data(), cell_specific_transformed_variables.data(), cell_specific_non_sampled_transformed_variables.data(), &value);
						return (Real)value;
					}
				}
			}
			return std::numeric_limits<Real>::quiet_NaN();
		}
	}
}

bool Cell::CellAliveAtTime(Real time, ESynchronizeCellTrajectory synchronize) const
{
	Real cell_time;
	if (synchronize == ESynchronizeCellTrajectory::DNAReplicationStart) {
		if (std::isnan(replication_start_time)) {
			cell_time = simulation_end_time;
		} else {
			cell_time = replication_start_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::PCNA_gfp_increase) {
		if (std::isnan(PCNA_gfp_increase_time)) {
			cell_time = simulation_end_time;
		} else {
			cell_time = PCNA_gfp_increase_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::NuclearEnvelopeBreakdown) {
		if (std::isnan(nuclear_envelope_breakdown_time)) {
			cell_time = simulation_end_time;
		} else {
			cell_time = nuclear_envelope_breakdown_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::AnaphaseOnset) {
		if (std::isnan(anaphase_onset_time)) {
			cell_time = simulation_end_time;
		} else {
			cell_time = anaphase_onset_time;
		}
	} else {
		cell_time = time - creation_time;
	}
	if (cell_time < 0.0 || cell_time > simulation_end_time) {
		return false;
	} else {
		return true;
	}
}

Real Cell::GetDuration(EPhaseDuration duration) const
{
	if (duration == EPhaseDuration::G1phase) {
		return replication_start_time;
	} else if (duration == EPhaseDuration::Sphase) {
		return replication_finish_time - replication_start_time;
	} else if (duration == EPhaseDuration::G2phase) {
		return nuclear_envelope_breakdown_time - replication_finish_time;
	} else if (duration == EPhaseDuration::NEBD_to_AnaphaseOnset) {
		return anaphase_onset_time - nuclear_envelope_breakdown_time;
	} else {
		LOGERROR("Unknown duration requested");
		return std::numeric_limits<Real>::quiet_NaN();
	}
}

void Cell::SetTreatmentConcentration(Real t)
{
	for (size_t i = 0; i < experiment->treatment_trajectories.size(); i++) {
		size_t ix = experiment->treatment_trajectories_species_ix[i];
		constant_species_y(ix) = experiment->treatment_trajectories[i]->GetConcentration(t, creation_time);
	}
}

bool Cell::solver_rhs_fn(OdeReal t, const OdeReal* y, OdeReal* ydot, void* user_data)
{
	SetTreatmentConcentration(t);
	if (Cell::use_generated_code) {
		experiment->derivative(ydot, y, constant_species_y.data(), cell_specific_transformed_variables.data(), cell_specific_non_sampled_transformed_variables.data());
	} else {
		model->CalculateDerivativePublic(t, y, ydot, constant_species_y.data(), cell_specific_transformed_variables.data(), cell_specific_non_sampled_transformed_variables.data());
		ASSERT(false);
	}
	return true;
}

bool Cell::solver_jac_fn(OdeReal t, const OdeReal* y, const OdeReal* ydot, OdeMatrixReal& jac, void* user_data)
{
	SetTreatmentConcentration(t);
	if (Cell::use_generated_code) {
		experiment->jacobian(jac, y, constant_species_y.data(), cell_specific_transformed_variables.data(), cell_specific_non_sampled_transformed_variables.data());
	} else {
		ASSERT(false);
	}
	return true;

}

Real Cell::discontinuity_cb(OdeReal t)
{
	Real disc = std::numeric_limits<Real>::infinity();
	for (size_t i = 0; i < experiment->treatment_trajectories.size(); i++) {
		Real d = experiment->treatment_trajectories[i]->NextDiscontinuity(t, creation_time);
		if (d < disc) {
			disc = d;
		}
	}
	if (disc == std::numeric_limits<Real>::infinity()) {
		return std::numeric_limits<Real>::quiet_NaN();
	} else {
		return disc;
	}
}

bool Cell::integration_step_cb(OdeReal t, const OdeReal* y, void* user_data)
{
	bool continue_integration = true;

	if (DNA_replication_ix != std::numeric_limits<size_t>::max() && replication_start_time != replication_start_time) {
		OdeReal s = y[DNA_replication_ix];
		if (s > 1e-4) {
			replication_start_time = solver->get_threshold_crossing_time(DNA_replication_ix, 1e-4, true, previous_integration_step_time);
		}
	}
	if (DNA_replicated_ix != std::numeric_limits<size_t>::max() && replication_finish_time != replication_finish_time) {
		OdeReal s = y[DNA_replicated_ix];
		if (s > 1.95) {
			replication_finish_time = solver->get_threshold_crossing_time(DNA_replicated_ix, 1.95, true, previous_integration_step_time);
		}
	}
	if (PCNA_gfp_ix != std::numeric_limits<size_t>::max() && PCNA_gfp_increase_time != PCNA_gfp_increase_time) {
		OdeReal s = y[PCNA_gfp_ix];
		if (s > 0.5) {
			PCNA_gfp_increase_time = solver->get_threshold_crossing_time(PCNA_gfp_ix, 0.5, true, previous_integration_step_time);
		}
	}
	if (nuclear_envelope_ix != std::numeric_limits<size_t>::max() && nuclear_envelope_breakdown_time != nuclear_envelope_breakdown_time) {
		OdeReal s = y[nuclear_envelope_ix];
		if (s < 0.5) {
			nuclear_envelope_breakdown_time = solver->get_threshold_crossing_time(nuclear_envelope_ix, 0.5, false, previous_integration_step_time);
		}
	}
	if (chromatid_separation_ix != std::numeric_limits<size_t>::max() && anaphase_onset_time != anaphase_onset_time) {
		if (y[chromatid_separation_ix] > 1e-3) {
			anaphase_onset_time = solver->get_threshold_crossing_time(chromatid_separation_ix, 1e-3, true, previous_integration_step_time);
			simulation_end_time = std::max(simulation_end_time, anaphase_onset_time + (OdeReal)simulate_past_chromatid_separation_time);
		}
	}
	if (cytokinesis_ix != std::numeric_limits<size_t>::max()) {
		if (y[cytokinesis_ix] > 1.0) {
			Real division_time = solver->get_threshold_crossing_time(cytokinesis_ix, 1.0, true, previous_integration_step_time);
			simulation_end_time = division_time;
			simulation_end_y = solver->GetInterpolatedY(division_time);
			if (experiment->divide_cells) {
				cell_divided = true;
				continue_integration = false;
			}
		}
	}
	if (apoptosis_ix != std::numeric_limits<size_t>::max()) {
		if (y[apoptosis_ix] > 1.0) {
			Real death_time = solver->get_threshold_crossing_time(apoptosis_ix, 1.0, true, previous_integration_step_time);
			simulation_end_time = death_time;
			simulation_end_y = solver->GetInterpolatedY(death_time);
			cell_died = true;
			continue_integration = false;
		}
	}

	previous_integration_step_time = t;

	return continue_integration;
}
