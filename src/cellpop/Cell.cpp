#include "Utils.h"
#include "Cell.h"
#include "ProbabilityDistributions.h"
#include "SBMLModel.h"
#include "SBMLSpecies.h"
#include "VariabilityDescription.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <fstream>
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"

size_t Cell::total_num_simulations = 0;
size_t Cell::cvode_max_steps_reached = 0;
size_t Cell::cvode_min_timestep_reached = 0;
const bool Cell::use_generated_code = 1;
static const int max_cvode_steps = 10000;

Cell::CVodeTimepoint::CVodeTimepoint()
	: cvode_time(0.0)
	, cv_uround(std::numeric_limits<Real>::quiet_NaN())
	, cv_tn(std::numeric_limits<Real>::quiet_NaN())
	, cv_h(std::numeric_limits<Real>::quiet_NaN())
	, cv_hu(std::numeric_limits<Real>::quiet_NaN())
	, cv_q(0)
{
}

static void static_cvode_err_fn(int error_code, const char *module, const char *function, char *msg, void *user_data)
{
//	LOGERROR("CVODE error %d in module %s, function %s: %s", error_code, module, function, msg);
}

int Cell::static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data)
{
	Cell* cell = reinterpret_cast<Cell*>(user_data);
	cell->SetTreatmentConcentration(t);
	if (Cell::use_generated_code) {
		cell->derivative(NV_DATA_S(ydot), NV_DATA_S(y), cell->constant_species_y.data(), cell->cell_specific_transformed_variables.data(), cell->cell_specific_non_sampled_transformed_variables.data());
	} else {
		cell->model->CalculateDerivativePublic(t, NV_DATA_S(y), NV_DATA_S(ydot), cell->constant_species_y.data(), cell->cell_specific_transformed_variables.data(), cell->cell_specific_non_sampled_transformed_variables.data());
	}

#if 0
	if (true) {
		for (int i = 0; i < NV_LENGTH_S(y); i++) {
			if (std::isnan(NV_Ith_S(y, i))) {
				int a;
				a = 6;
			}
			if (std::isnan(NV_Ith_S(ydot, i))) {
				int a;
				a = 6;
			}
		}
		std::ofstream file("tmp.txt", std::ios::app);
		file.precision(12);
		file << t << "\t";
		for (int i = 0; i < NV_LENGTH_S(y); i++) {
			file << NV_Ith_S(y, i) << "\t";
		}
		for (int i = 0; i < NV_LENGTH_S(ydot); i++) {
			file << NV_Ith_S(ydot, i) << "\t";
		}
		file << std::endl;
	}
#endif

	return 0;
}

int Cell::static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	Cell* cell = reinterpret_cast<Cell*>(user_data);
	cell->SetTreatmentConcentration(t);
	cell->jacobian(Jac, NV_DATA_S(y), cell->constant_species_y.data(), cell->cell_specific_transformed_variables.data(), cell->cell_specific_non_sampled_transformed_variables.data());

#if 0
	for (int i = 0; i < NV_LENGTH_S(y); i++) {
		for (int j = 0; j < NV_LENGTH_S(y); j++) {
			if (std::isnan(SM_ELEMENT_D(Jac, i, j))) {
				int a;
				a = 6;
			}
		}
	}

	std::ofstream file("tmp_jac.txt", std::ios::app);
	file.precision(18);
	file << "Jacobian; t=" << t << std::endl;
	file << "y=" << EIGV(y).transpose() << std::endl;
	file << EIGMAT(Jac) << std::endl;
#endif

	return 0;
}

Cell::Cell(const SBMLModel* model, const Experiment* experiment)
	: model(model)
	, experiment(experiment)
	, cvode_initialized(false)
	, creation_time(std::numeric_limits<Real>::quiet_NaN())
	, current_simulation_time(std::numeric_limits<Real>::quiet_NaN())
	, cvode_steps(0)
	, min_step_size(std::numeric_limits<Real>::infinity())
	, completed(false)
	, interpolation_time(-std::numeric_limits<Real>::infinity())
	, derivative(NULL)
	, jacobian(NULL)
	, synchronize_offset_time(0.0)
	, cvode_timepoint_iter(0)
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
	size_t num_cvode_species = model->GetNumCVodeSpecies();

	cvode_mem = CVodeCreate(CV_BDF);
	cvode_y = MakeCVodeVector(num_cvode_species);
	cvode_interpolate_y.setConstant(num_cvode_species, 0.0);
	J = MakeCVodeMatrix(num_cvode_species, num_cvode_species);
	LS = MakeCVodeLinearSolver(cvode_y, J);
	NLS = SUNNonlinSol_Newton(cvode_y);

	constant_species_y.setZero(model->GetNumSimulatedSpecies() - num_cvode_species);
	cell_specific_transformed_variables.setConstant(experiment->GetVarset()->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());
	cell_specific_non_sampled_transformed_variables.setConstant(experiment->non_sampled_parameters.size(), std::numeric_limits<Real>::quiet_NaN());

	cvode_timepoints.resize(max_cvode_steps);
	for (size_t j = 0; j < 6; j++) {
		cvode_timepoints_zn[j].resize(num_cvode_species, max_cvode_steps);
		//cvode_timepoints_zn[j].setConstant(num_cvode_species, max_cvode_steps, std::numeric_limits<Real>::quiet_NaN());
	}

	DNA_replication_ix = model->GetCVodeSpeciesByName("replicating_DNA", false);
	DNA_replicated_ix = model->GetCVodeSpeciesByName("replicated_DNA", false);
	PCNA_gfp_ix = model->GetCVodeSpeciesByName("PCNA_gfp", false);
	nuclear_envelope_ix = model->GetCVodeSpeciesByName("nuclear_envelope", false);
	chromatid_separation_ix = model->GetCVodeSpeciesByName("chromatid_separation", false);
	cytokinesis_ix = model->GetCVodeSpeciesByName("cytokinesis", false);
	apoptosis_ix = model->GetCVodeSpeciesByName("apoptosis", false);
}

Cell::~Cell()
{
	SUNMatDestroy(J);
	N_VDestroy(cvode_y);
	SUNNonlinSolFree(NLS);
	SUNLinSolFree(LS);
	CVodeFree(&cvode_mem);
}

void Cell::SetDerivativeFunctions(derivative_fn fn, jacobian_fn jac)
{
	derivative = fn;
	jacobian = jac;
}

bool Cell::SetInitialConditionsFromModel(const std::map<size_t, Experiment::SetSpecies>& set_species_map, const std::map<size_t, size_t>& set_init_map, const std::map<size_t, std::vector<int>>& ratio_active_map,const std::map<size_t, std::vector<int>>& ratio_inactive_map, const std::map<size_t, std::vector<size_t>>& ratio_total_active, const std::map<size_t, std::vector<size_t>>& ratio_total_inactive,const VectorReal& transformed_values, Real time)
{
	for (size_t i = 0; i < model->GetNumCVodeSpecies(); i++) {
		NV_Ith_S(cvode_y, i) = model->GetCVodeSpecies(i)->GetInitialValue();
	}
	for (size_t i = 0; i < model->GetNumConstantSpecies(); i++) {
		constant_species_y(i) = model->GetConstantSpecies(i)->GetInitialValue();
	}

	for (std::map<size_t, Experiment::SetSpecies>::const_iterator ssmi = set_species_map.begin(); ssmi != set_species_map.end(); ++ssmi) {
		if (time >= ssmi->second.begin_time && time < ssmi->second.end_time) {
			NV_Ith_S(cvode_y, ssmi->first) = ssmi->second.value;
		}
	}

	for(auto const& sic : set_init_map){
		NV_Ith_S(cvode_y, sic.first) = transformed_values[sic.second];
	}

	for(auto const& ram : ratio_active_map){
		NV_Ith_S(cvode_y, ram.first) = transformed_values[ram.second[0]] * transformed_values[ram.second[1]];
	}

	for(auto const& rim : ratio_inactive_map){
		NV_Ith_S(cvode_y, rim.first) = (1 - transformed_values[rim.second[0]]) * transformed_values[rim.second[1]];
	}

	for(auto const& rtm : ratio_total_active){
		NV_Ith_S(cvode_y, rtm.first) = transformed_values[rtm.second[0]] * (model -> GetCVodeSpecies(rtm.second[1]) -> GetInitialValue() + model -> GetCVodeSpecies(rtm.second[2]) -> GetInitialValue());
	}

	for(auto const& rtm : ratio_total_inactive){
		NV_Ith_S(cvode_y, rtm.first) = (1 - transformed_values[rtm.second[0]]) * (model -> GetCVodeSpecies(rtm.second[1]) -> GetInitialValue() + model -> GetCVodeSpecies(rtm.second[2]) -> GetInitialValue());
	}

	return true;
}

bool Cell::SetInitialConditionsFromOtherCell(const Cell* other)
{
	for (size_t i = 0; i < model->GetNumCVodeSpecies(); i++) {
		NV_Ith_S(cvode_y, i) = other->cvode_end_y(i);
	}
	constant_species_y = other->constant_species_y;

#if 1
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("DNA", true)) *= 0.5;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("replicating_DNA", true)) *= 0.5;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("replicated_DNA", true)) *= 0.5;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("licensed_DNA", true)) *= 0.5;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("mitogenic_signal", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("unassembled_spindle", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("assembled_spindle", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("chromatid_separation", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("nuclear_envelope", true)) = 1.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("cytokinesis", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("G2_delay", true)) = 0.0;
#else
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("cell_cycle_start", true)) = 1.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("G1_checkpoint", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("S-G2-M_phase", true)) = 0.0;
	NV_Ith_S(cvode_y, model->GetCVodeSpeciesByName("cytokinesis", true)) = 0.0;
#endif

	return true;
}

bool Cell::Initialize(Real creation_time, const VectorReal& transformed_variables, VectorReal* sobol_sequence_values, bool is_initial_cell, bool calculate_synchronization_point, Real abs_tol, Real rel_tol)
{
	cvode_steps = 0;
	cvode_timepoint_iter = 0;

	int sobol_sequence_ix = 0;
	for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
		(*it)->ApplyVariabilityEntryTime(creation_time, *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
	}
	this->creation_time = creation_time;

	for (size_t i = 0; i < cell_specific_transformed_variables.size(); i++) {
		cell_specific_transformed_variables(i) = transformed_variables(i);
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilityParameter(experiment->varset->GetVariableName(i), cell_specific_transformed_variables(i), *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
		}
	}

	for (size_t i = 0; i < cell_specific_non_sampled_transformed_variables.size(); i++) {
		cell_specific_non_sampled_transformed_variables(i) = experiment->non_sampled_parameters(i);
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilityParameter(experiment->non_sampled_parameter_names[i], cell_specific_non_sampled_transformed_variables(i), *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
		}
	}

	for (size_t i = 0; i < model->GetNumCVodeSpecies(); i++) {
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilitySpecies(model->GetCVodeSpeciesName(i), NV_Ith_S(cvode_y, i), *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters, is_initial_cell);
		}
	}

	completed = false;
	replication_start_time = std::numeric_limits<Real>::quiet_NaN();
	replication_finish_time = std::numeric_limits<Real>::quiet_NaN();
	PCNA_gfp_increase_time = std::numeric_limits<Real>::quiet_NaN();
	nuclear_envelope_breakdown_time = std::numeric_limits<Real>::quiet_NaN();
	anaphase_onset_time = std::numeric_limits<Real>::quiet_NaN();

	if (!cvode_initialized) {
		CVodeInit(cvode_mem, static_cvode_rhs_fn, 0.0, cvode_y);
		cvode_initialized = true;
	} else {
		CVodeReInit(cvode_mem, 0.0, cvode_y);
	}
	CVodeSetUserData(cvode_mem, (void*)this);
	CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);
	CVodeSStolerances(cvode_mem, abs_tol, rel_tol);
	CVodeSetMinStep(cvode_mem, 1e-3);
	//CVodeSetMaxStep(cvode_mem, 30.0*60.0);
	CVodeSetMaxStep(cvode_mem, 5.0);
	CVodeSetInitStep(cvode_mem, 1.0);

	CVodeSetLinearSolver(cvode_mem, LS, J);

	//FILE* infofp = fopen("cvode_monitor.txt", "w+");
	//SUNNonlinSolSetPrintLevel_Newton(NLS, 1);
	//SUNNonlinSolSetInfoFile_Newton(NLS, infofp);
	CVodeSetNonlinearSolver(cvode_mem, NLS);

	if (Cell::use_generated_code) {
		CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);
		CVodeSetMaxStepsBetweenJac(cvode_mem, 10);
	}

	total_num_simulations++;

	return true;
}

bool Cell::Simulate(Real end_time, Real simulate_past_chromatid_separation_time, bool& die, bool& divide, Real& achieved_time)
{
	die = false;
	divide = false;
	achieved_time = 0.0;

	current_simulation_time = 0;
	bool result = true;
	OdeReal cell_end_time = (OdeReal)end_time - creation_time;
	while (current_simulation_time < cell_end_time) {
		if (cvode_steps >= max_cvode_steps) {
			cvode_max_steps_reached++;
			//printf("CVode max steps");
			result = false;
			break;
		}

		OdeReal prev_time = current_simulation_time;
		int cvode_result = CVode(cvode_mem, cell_end_time, cvode_y, &current_simulation_time, CV_ONE_STEP);
		if (cvode_result < 0) {
			// Try once more
			cvode_result = CVode(cvode_mem, cell_end_time, cvode_y, &current_simulation_time, CV_ONE_STEP);
		}
		if (cvode_result == CV_ERR_FAILURE || cvode_result == CV_CONV_FAILURE) {
			cvode_min_timestep_reached++;
		}
		if (cvode_result < 0) {
			//printf("CVode failure: %u", cvode_result);
			result = false;
			break;
		}

#if 0
		OdeReal hcur;
		CVodeGetCurrentStep(cvode_mem, &hcur);
		int qcur;
		CVodeGetCurrentOrder(cvode_mem, &qcur);

		std::ofstream file("cvode_stats.txt", std::ios::app);
		file.precision(18);
		file << "CVode step " << std::to_string(cvode_steps + 1);
		file << "; time=" << current_simulation_time;
		file << " hcur=" << hcur;
		file << " qcur=" << qcur;
		file <<std::endl;
		//file << "cvode_y=" << EIGV(cvode_y).transpose() << std::endl;
		file.close();

#if CVODE_USE_EIGEN_SOLVER
		{
			std::ofstream file("trajectory.txt", std::ios::app);
			file.precision(18);
			file << "" << current_simulation_time;
			file << ";" << EIGV(cvode_y).transpose();
			OdeVectorReal& acor = EIGV(((CVodeMem)cvode_mem)->cv_acor);
			file << ";" << acor.transpose();
			file << std::endl;
			//file << "cvode_y=" << EIGV(cvode_y).transpose() << std::endl;
			file.close();
		}
#endif
#endif

		// Store relevant information for interpolation at any timepoint later
		RetrieveCVodeInterpolationInfo();
		min_step_size = (std::min)(min_step_size, current_simulation_time - prev_time);

		if (DNA_replication_ix != std::numeric_limits<size_t>::max() && replication_start_time != replication_start_time) {
			OdeReal s = NV_Ith_S(cvode_y, DNA_replication_ix);
			if (s > 1e-4) {
				replication_start_time = InterpolateEventTime(DNA_replication_ix, 1e-4, true, prev_time);
			}
		}
		if (DNA_replicated_ix != std::numeric_limits<size_t>::max() && replication_finish_time != replication_finish_time) {
			OdeReal s = NV_Ith_S(cvode_y, DNA_replicated_ix);
			if (s > 1.95) {
				replication_finish_time = InterpolateEventTime(DNA_replicated_ix, 1.95, true, prev_time);
			}
		}
		if (PCNA_gfp_ix != std::numeric_limits<size_t>::max() && PCNA_gfp_increase_time != PCNA_gfp_increase_time) {
			OdeReal s = NV_Ith_S(cvode_y, PCNA_gfp_ix);
			if (s > 0.5) {
				PCNA_gfp_increase_time = InterpolateEventTime(PCNA_gfp_ix, 0.5, true, prev_time);
			}
		}
		if (nuclear_envelope_ix != std::numeric_limits<size_t>::max() && nuclear_envelope_breakdown_time != nuclear_envelope_breakdown_time) {
			OdeReal s = NV_Ith_S(cvode_y, nuclear_envelope_ix);
			if (s < 0.5) {
				nuclear_envelope_breakdown_time = InterpolateEventTime(nuclear_envelope_ix, 0.5, false, prev_time);
			}
		}
		if (chromatid_separation_ix != std::numeric_limits<size_t>::max() && anaphase_onset_time != anaphase_onset_time) {
			if (NV_Ith_S(cvode_y, chromatid_separation_ix) > 1e-3) {
				anaphase_onset_time = InterpolateEventTime(chromatid_separation_ix, 1e-3, true, prev_time);
				cell_end_time = std::max(cell_end_time, anaphase_onset_time + (OdeReal)simulate_past_chromatid_separation_time);
			}
		}
		if (cytokinesis_ix != std::numeric_limits<size_t>::max()) {
			if (NV_Ith_S(cvode_y, cytokinesis_ix) > 1.0) {
				achieved_time = InterpolateEventTime(cytokinesis_ix, 1.0, true, prev_time);
				CalculateEndY(achieved_time);
				divide = true;
			}
		}
		if (apoptosis_ix != std::numeric_limits<size_t>::max()) {
			if (NV_Ith_S(cvode_y, apoptosis_ix) > 1.0) {
				achieved_time = InterpolateEventTime(apoptosis_ix, 1.0, true, prev_time);
				CalculateEndY(achieved_time);
				die = true;
			}
		}

		if (die || (experiment->divide_cells && divide)) {
			completed = true;
			break;
		}

		cvode_steps++;
	}

	//printf("Simulation steps: %zu\n", cvode_steps);
	achieved_time = current_simulation_time + creation_time;
	return result;
}

Real Cell::GetInterpolatedSpeciesValue(Real time, size_t species_ix, ESynchronizeCellTrajectory synchronize)
{
	if (species_ix >= model->GetNumCVodeSpecies()) {
		LOGERROR("Out of bounds species index");
		ASSERT(false);
		return std::numeric_limits<Real>::quiet_NaN();
	}
	if (cvode_steps == 0) {
		return std::numeric_limits<Real>::quiet_NaN();
	}

	Real cell_time = 0.0;
	if (synchronize == ESynchronizeCellTrajectory::DNAReplicationStart) {
		if (std::isnan(replication_start_time)) {
			cell_time = time + cvode_timepoints[cvode_steps - 1].cvode_time;
		} else {
			cell_time = time + replication_start_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::PCNA_gfp_increase) {
		if (std::isnan(PCNA_gfp_increase_time)) {
			cell_time = time + cvode_timepoints[cvode_steps - 1].cvode_time;
		} else {
			cell_time = time + PCNA_gfp_increase_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::NuclearEnvelopeBreakdown) {
		if (std::isnan(nuclear_envelope_breakdown_time)) {
			cell_time = time + cvode_timepoints[cvode_steps-1].cvode_time;
		} else {
			cell_time = time + nuclear_envelope_breakdown_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::AnaphaseOnset) {
		if (std::isnan(anaphase_onset_time)) {
			cell_time = time + cvode_timepoints[cvode_steps-1].cvode_time;
		} else {
			cell_time = time + anaphase_onset_time;
		}
	} else {
		cell_time = time - creation_time;
	}
	if (cell_time < 0.0) {
		// This cell didn't exist yet..
		return std::numeric_limits<Real>::quiet_NaN();
	}

	if (cell_time == interpolation_time) {
		// We still have this interpolation
		return cvode_interpolate_y(species_ix);
	}

	// Requests for interpolations should be in strictly increasing time, so we can continue searching where we left off
	while (cvode_timepoint_iter < cvode_steps) {
		const CVodeTimepoint& cvt = cvode_timepoints[cvode_timepoint_iter];
		if (cvt.cvode_time > cell_time) {
			break;
		} else {
			cvode_timepoint_iter++;
		}
	}
	if (cvode_timepoint_iter == cvode_steps) {
		// This timepoint is beyond the cell's simulation time; apparently the cell died or divided
		return std::numeric_limits<Real>::quiet_NaN();
	}
	interpolation_time = cell_time;

	// Calculate the interpolated value
	// This is based on CVodeGetDky with k=0; we've stored the relevant part of cvode_mem in the CVodeTimepoint struct
	const CVodeTimepoint& cvt = cvode_timepoints[cvode_timepoint_iter];

	/* Allow for some slack */
	Real tfuzz = 100.0 * cvt.cv_uround * fabs(cvt.cv_tn) + fabs(cvt.cv_hu);
	if (cvt.cv_hu < 0.0) {
		tfuzz = -tfuzz;
	}
	Real tp = cvt.cv_tn - cvt.cv_hu - tfuzz;
	Real tn1 = cvt.cv_tn + tfuzz;
	if ((cell_time - tp) * (cell_time - tn1) > 0.0) {
		LOGERROR("Time error for interpolation");
		return std::numeric_limits<Real>::quiet_NaN();
	}

	/* Sum the differentiated interpolating polynomial */
	CVodeMem cv_mem = (CVodeMem)cvode_mem;
	Real s = (cell_time - cvt.cv_tn) / cvt.cv_h;
	cvode_interpolate_y.setZero();
	for (int j = cvt.cv_q; j >= 0; j--) {
		Real cval = 1.0;
		for (int i = j; i >= j + 1; i--)
			cval *= i;
		for (int i = 0; i < j; i++)
			cval *= s;

		cvode_interpolate_y += cval * cvode_timepoints_zn[j].col(cvode_timepoint_iter);
	}

	return cvode_interpolate_y(species_ix);
}

void Cell::RestartInterpolationIteration()
{
	cvode_interpolate_y.setZero();
	cvode_timepoint_iter = 0;
	interpolation_time = std::numeric_limits<Real>::quiet_NaN();
}

bool Cell::CellAliveAtTime(Real time, ESynchronizeCellTrajectory synchronize) const
{
	if (cvode_steps == 0) {
		return false;
	}
	Real cell_time = 0.0;
	if (synchronize == ESynchronizeCellTrajectory::DNAReplicationStart) {
		if (std::isnan(replication_start_time)) {
			cell_time += cvode_timepoints[cvode_steps - 1].cvode_time;
		} else {
			cell_time += replication_start_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::PCNA_gfp_increase) {
		if (std::isnan(PCNA_gfp_increase_time)) {
			cell_time += cvode_timepoints[cvode_steps - 1].cvode_time;
		} else {
			cell_time += PCNA_gfp_increase_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::NuclearEnvelopeBreakdown) {
		if (std::isnan(nuclear_envelope_breakdown_time)) {
			cell_time += cvode_timepoints[cvode_steps - 1].cvode_time;
		} else {
			cell_time += nuclear_envelope_breakdown_time;
		}
	} else if (synchronize == ESynchronizeCellTrajectory::AnaphaseOnset) {
		if (std::isnan(anaphase_onset_time)) {
			cell_time += cvode_timepoints[cvode_steps - 1].cvode_time;
		} else {
			cell_time += anaphase_onset_time;
		}
	} else {
		cell_time = time - creation_time;
	}
	if (cell_time < 0.0) {
		return false;
	}
	if (cell_time > cvode_timepoints[cvode_steps-1].cvode_time) {
		return false;
	}
	return true;
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

void Cell::SetMutations()
{
	// hack
	//constant_species_y(2) = 1.0;
}

inline Real pulse(Real t, Real t_in)
{
	Real t_in_pulse = t - t_in - 2.0;

	//return ((0.5 * tanh(100.0 * t_in - 100.0 * t + 300.0) + 0.5) * (t_in - 1.0 * t + 2.0) + (0.5 * tanh(100.0 * t_in - 100.0 * t + 300.0) - 0.5) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1000.0) - 1.0 * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1400.0) + 0.5) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1000.0) - 0.5) * (0.25 * t_in - 0.25 * t + 3.5) + 0.5)) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 200.0) - 0.5);
	if (t_in_pulse < 0.0) {
		return 0.0;
	} else if (t_in_pulse < 2.0) {
		return t_in_pulse * 0.5;
	} else if (t_in_pulse < 10.0) {
		return 1.0;
	} else if (t_in_pulse < 14.0) {
		return 1 - (t_in_pulse - 10.0) * 0.25;
	} else {
		return 0.0;
	}
}

void Cell::SetTreatmentConcentration(Real t)
{
	Real t_in = 30.0;
	//Real test = ((0.5 * tanh(100.0 * t_in - 100.0 * t + 300.0) + 0.5) * (t_in - 1.0 * t + 2.0) + (0.5 * tanh(100.0 * t_in - 100.0 * t + 300.0) - 0.5) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1000.0) - 1.0 * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1400.0) + 0.5) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1000.0) - 0.5) * (0.25 * t_in - 0.25 * t + 3.5) + 0.5)) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 200.0) - 0.5);
	constant_species_y(0) = 0.0;
	if (t < 60) constant_species_y(0) += pulse(t, 30);
	else if (t < 90) constant_species_y(0) += pulse(t, 60);
	else if (t < 120) constant_species_y(0) += pulse(t, 90);
	else if (t < 150) constant_species_y(0) += pulse(t, 120);
	else if (t < 180) constant_species_y(0) += pulse(t, 150);
	else if (t < 210) constant_species_y(0) += pulse(t, 180);
	else if (t < 241.9064) constant_species_y(0) += pulse(t, 210);
	else if (t < 290.4352) constant_species_y(0) += pulse(t, 241.9064);
	else if (t < 372.1597) constant_species_y(0) += pulse(t, 290.4352);
	else if (t < 455.319) constant_species_y(0) += pulse(t, 372.1597);
	else if (t < 510.8225) constant_species_y(0) += pulse(t, 455.319);
	else if (t < 546.0478) constant_species_y(0) += pulse(t, 510.8225);
	else constant_species_y(0) += pulse(t, 546.0478);

	for (size_t i = 0; i < experiment->treatment_trajectories.size(); i++) {
		size_t ix = experiment->treatment_trajectories_species_ix[i];
		constant_species_y(ix) = experiment->treatment_trajectories[i]->GetConcentration(t + creation_time, experiment->selected_treatment_trajectory_sample[i]);
	}
}

void Cell::RetrieveCVodeInterpolationInfo()
{
	CVodeTimepoint& cvt = cvode_timepoints[cvode_steps];

	// Copy relevant variables from cvode_mem
	CVodeMem cv_mem = (CVodeMem)cvode_mem;

	cvt.cvode_time	= current_simulation_time;
	cvt.cv_uround	= cv_mem->cv_uround;
	cvt.cv_tn		= cv_mem->cv_tn;
	cvt.cv_h		= cv_mem->cv_h;
	cvt.cv_hu		= cv_mem->cv_hu;
	cvt.cv_q		= cv_mem->cv_q;

	ASSERT(cvt.cv_q <= 5);
	for (int j = cv_mem->cv_q; j >= 0; j--) {
#if CVODE_USE_EIGEN_SOLVER
		cvode_timepoints_zn[j].col(cvode_steps) = *NV_CONTENT_S(cv_mem->cv_zn[j]);
#else
		for (int i = 0; i < NV_LENGTH_S(cv_mem->cv_zn[j]); i++) {
			cvode_timepoints_zn[j](i, cvode_steps) = NV_Ith_S(cv_mem->cv_zn[j], i);
		}
#endif
	}
}

Real Cell::InterpolateEventTime(size_t species_ix, Real threshold, bool above, Real prev_time)
{
	CVodeTimepoint& cvt = cvode_timepoints[cvode_steps];

	Real dt = (current_simulation_time - prev_time) * 0.5;
	Real time = prev_time + dt;

	bool crossed = true;
	for (int iter = 0; iter < 10; iter++) {
		/* Allow for some slack */
		Real tfuzz = 100.0 * cvt.cv_uround * fabs(cvt.cv_tn) + fabs(cvt.cv_hu);
		if (cvt.cv_hu < 0.0) {
			tfuzz = -tfuzz;
		}
		Real tp = cvt.cv_tn - cvt.cv_hu - tfuzz;
		Real tn1 = cvt.cv_tn + tfuzz;

		/* Sum the differentiated interpolating polynomial */
		CVodeMem cv_mem = (CVodeMem)cvode_mem;
		Real s = (time - cvt.cv_tn) / cvt.cv_h;
		Real x = 0.0;
		for (int j = cvt.cv_q; j >= 0; j--) {
			Real cval = 1.0;
			for (int i = j; i >= j + 1; i--)
				cval *= i;
			for (int i = 0; i < j; i++)
				cval *= s;

			x += cval * cvode_timepoints_zn[j](species_ix, cvode_steps);
		}

		dt *= 0.5;
		if (above) {
			if (x > threshold) {
				time -= dt;
			} else {
				time += dt;
			}
		} else {
			if (x < threshold) {
				time -= dt;
			} else {
				time += dt;
			}
		}
	}
	
	return time;
}

void Cell::CalculateEndY(Real end_time)
{
	CVodeTimepoint& cvt = cvode_timepoints[cvode_steps];

	cvode_end_y.setConstant(model->GetNumCVodeSpecies(), 0.0);

	/* Allow for some slack */
	Real tfuzz = 100.0 * cvt.cv_uround * fabs(cvt.cv_tn) + fabs(cvt.cv_hu);
	if (cvt.cv_hu < 0.0) {
		tfuzz = -tfuzz;
	}
	Real tp = cvt.cv_tn - cvt.cv_hu - tfuzz;
	Real tn1 = cvt.cv_tn + tfuzz;

	/* Sum the differentiated interpolating polynomial */
	CVodeMem cv_mem = (CVodeMem)cvode_mem;
	Real s = (end_time - cvt.cv_tn) / cvt.cv_h;
	Real x = 0.0;
	for (int j = cvt.cv_q; j >= 0; j--) {
		Real cval = 1.0;
		for (int i = j; i >= j + 1; i--)
			cval *= i;
		for (int i = 0; i < j; i++)
			cval *= s;

		cvode_end_y += cval * cvode_timepoints_zn[j].col(cvode_steps);
	}
}
