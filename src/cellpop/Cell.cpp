#include "Utils.h"
#include "Cell.h"
#include "ProbabilityDistributions.h"
#include "SBMLModel.h"
#include "SBMLSpecies.h"
#include "VariabilityDescription.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"

#define USE_GENERATED_CODE 1

size_t Cell::total_num_simulations = 0;
size_t Cell::cvode_max_steps_reached = 0;
size_t Cell::cvode_min_timestep_reached = 0;
static const int max_cvode_steps = 5000;

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

bool tmp = false;
int Cell::static_cvode_rhs_fn(Real t, N_Vector y, N_Vector ydot, void* user_data)
{
	Cell* cell = reinterpret_cast<Cell*>(user_data);
	cell->SetTreatmentConcentration(t);
#if USE_GENERATED_CODE
	cell->derivative(NV_DATA_S(ydot), NV_DATA_S(y), cell->constant_species_y.data(), cell->cell_specific_transformed_variables.data(), cell->cell_specific_non_sampled_transformed_variables.data());
#else
	cell->model->CalculateDerivativePublic(t, NV_DATA_S(y), NV_DATA_S(ydot), cell->constant_species_y.data(), cell->cell_specific_transformed_variables.data());
#endif

#if 0
	if (t == 0.0 && !tmp) {
		for (int i = 0; i < NV_LENGTH_S(y); i++) {
			LOG("u0[%d] = %.20g;", i, NV_Ith_S(y, i));
		}
		for (int i = 0; i < cell->cell_specific_transformed_variables.size(); i++) {
			LOG("params[%d] = %.20g;", i, cell->cell_specific_transformed_variables(i));
		}
		for (int i = 0; i < cell->constant_species_y.size(); i++) {
			LOG("constant_species[%d] = %.20g;", i, cell->constant_species_y(i));
		}
		for (int i = 0; i < NV_LENGTH_S(y); i++) {
			LOG("du[%d] = %.20g;", i, NV_Ith_S(ydot, i));
		}
		tmp = true;
	}
	for (size_t i = 0; i < NV_LENGTH_S(ydot); i++) {
		if (std::isnan(NV_Ith_S(ydot, i))) {
			int a;
			a = 9;
		}
		if (NV_Ith_S(ydot, i) != 0.0 && !std::isnormal(NV_Ith_S(ydot, i))) {
			int a;
			a = 9;
		}
		if (NV_Ith_S(y, i) != 0.0 && !std::isnormal(NV_Ith_S(y, i))) {
			int a;
			a = 9;
		}
		if (NV_Ith_S(y, i) <= -1e-4 && NV_Ith_S(ydot, i) < 0.0) {
			int a;
			a = 9;
		}
	}
#endif

	return 0;
}

#if USE_GENERATED_CODE
static int tmpvar = 0;
bool tmp2 = false;
int Cell::static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	Cell* cell = reinterpret_cast<Cell*>(user_data);
	cell->SetTreatmentConcentration(t);
	//return solver->cvode_jac_fn(t, NV_DATA_S(y), NV_DATA_S(fy), Jac->cols);
	cell->jacobian(Jac, NV_DATA_S(y), cell->constant_species_y.data(), cell->cell_specific_transformed_variables.data(), cell->cell_specific_non_sampled_transformed_variables.data());

#if 0
	if (!tmp2) {
		for (int i = 0; i < NV_LENGTH_S(y); i++) {
			LOG("u0[%d] = %.20g;", i, NV_Ith_S(y, i));
		}
		for (int i = 0; i < cell->cell_specific_transformed_variables.size(); i++) {
			LOG("params[%d] = %.20g;", i, cell->cell_specific_transformed_variables(i));
		}
		for (int i = 0; i < cell->constant_species_y.size(); i++) {
			LOG("constant_species[%d] = %.20g;", i, cell->constant_species_y(i));
		}
		for (int i = 0; i < 13; i++) {
			for (int j = 0; j < 13; j++) {
				if (SM_ELEMENT_D(Jac, i, j) != 0.0) {
					LOG("out[%d,%d] = %.20g;", i, j, SM_ELEMENT_D(Jac, i, j));
				}
			}
		}
		tmp2 = true;
	}
	//if (tmpvar == 0) {
		printf("%d\n", tmpvar);
		for (size_t i = 0; i < NV_LENGTH_S(y); i++) {
			for (size_t j = 0; j < NV_LENGTH_S(y); j++) {
				printf("%12g,", SM_ELEMENT_D(Jac, i, j));
				//			if (SM_ELEMENT_D(Jac, i, j) != SM_ELEMENT_D(Jac, i, j)) {
				//				int a;
				//				a = 9;
				//			}
			}
			printf("\n");
		}
		tmpvar++;
	//}
#endif

	return 0;
}
#endif

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
#if 0
	long int nsteps;
	long int nfevals;
	long int nlinsetups;
	long int netfails;
	int qlast;
	int qcur;
	realtype hinused;
	realtype hlast;
	realtype hcur;
	realtype tcur;
	CVodeGetIntegratorStats(cvode_mem, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast, &qcur, &hinused, &hlast, &hcur, &tcur);
	printf("%d %d %d %d %d %d %g %g %g %g", nsteps, nfevals, nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur);
#endif

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

bool Cell::SetInitialConditionsFromModel(const std::map<size_t, Experiment::SetSpecies>& set_species_map, Real time)
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

	return true;
}

bool Cell::SetInitialConditionsFromOtherCell(const Cell* other)
{
	// There's no direct set function?
	N_VAddConst(other->cvode_y, 0.0, cvode_y);
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

bool Cell::Initialize(Real creation_time, const VectorReal& transformed_variables, VectorReal* sobol_sequence_values, bool apply_entry_time_variability, bool calculate_synchronization_point)
{
	cvode_steps = 0;
	cvode_timepoint_iter = 0;

	int sobol_sequence_ix = 0;
	if (apply_entry_time_variability) {
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilityEntryTime(creation_time, *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters);
		}
	} else {
		// Discard as many random variables from the sobol sequence as necessary to account for the fact
		// that this cell does not have variable entry time
		Real dummy = 0.0;
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilityEntryTime(dummy, *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters);
		}
	}
	this->creation_time = creation_time;

	for (size_t i = 0; i < cell_specific_transformed_variables.size(); i++) {
		cell_specific_transformed_variables(i) = transformed_variables(i);
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilityParameter(experiment->varset->GetVariableName(i), cell_specific_transformed_variables(i), *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters);
		}
	}

	for (size_t i = 0; i < cell_specific_non_sampled_transformed_variables.size(); i++) {
		cell_specific_non_sampled_transformed_variables(i) = experiment->non_sampled_parameters(i);
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilityParameter(experiment->non_sampled_parameter_names[i], cell_specific_non_sampled_transformed_variables(i), *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters);
		}
	}

	for (size_t i = 0; i < model->GetNumCVodeSpecies(); i++) {
		for (auto it = experiment->cell_variabilities.begin(); it != experiment->cell_variabilities.end(); ++it) {
			(*it)->ApplyVariabilitySpecies(model->GetCVodeSpeciesName(i), NV_Ith_S(cvode_y, i), *sobol_sequence_values, sobol_sequence_ix, transformed_variables, experiment->non_sampled_parameters);
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
	CVodeSStolerances(cvode_mem, 1e-8, 1e-8);
	//CVodeSetMinStep(cvode_mem, 1e-3);
	CVodeSetInitStep(cvode_mem, 1.0);

	CVodeSetLinearSolver(cvode_mem, LS, J);

	//FILE* infofp = fopen("cvode_monitor.txt", "w+");
	//SUNNonlinSolSetPrintLevel_Newton(NLS, 1);
	//SUNNonlinSolSetInfoFile_Newton(NLS, infofp);
	CVodeSetNonlinearSolver(cvode_mem, NLS);

#if USE_GENERATED_CODE
	CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);
	CVodeSetMaxStepsBetweenJac(cvode_mem, 10);
#endif

	total_num_simulations++;

	return true;
}

#if 0
bool Cell::SetSpecies(const std::map<size_t, Real>& set_species_map)
{
	for (std::map<size_t, Real>::const_iterator ssmi = set_species_map.begin(); ssmi != set_species_map.end(); ++ssmi) {
		NV_Ith_S(cvode_y, ssmi->first) = ssmi->second;
	}
	CVodeReInit(cvode_mem, current_simulation_time, cvode_y);
	return true;
}

bool Cell::SetSpecies(size_t i, Real value, bool reinit)
{
	NV_Ith_S(cvode_y, i) = value;
	if (reinit) {
		CVodeReInit(cvode_mem, current_simulation_time, cvode_y);
	}
	return true;
}
#endif

bool Cell::Simulate(Real end_time, bool& die, bool& divide, Real& achieved_time)
{
	die = false;
	divide = false;
	achieved_time = 0.0;

	current_simulation_time = 0;
	bool result = true;
	Real cell_end_time = end_time - creation_time;
	while (current_simulation_time < cell_end_time) {
		if (cvode_steps >= max_cvode_steps) {
			cvode_max_steps_reached++;
			//printf("CVode max steps");
			result = false;
			break;
		}

		Real prev_time = current_simulation_time;
		int result = CVode(cvode_mem, cell_end_time, cvode_y, &current_simulation_time, CV_ONE_STEP);
		if (result < 0) {
			// Try once more
			result = CVode(cvode_mem, cell_end_time, cvode_y, &current_simulation_time, CV_ONE_STEP);
		}
		if (result == CV_ERR_FAILURE || result == CV_CONV_FAILURE) {
			cvode_min_timestep_reached++;
		}
		if (result < 0) {
			//printf("CVode failure: %u", result);
			result = false;
			break;
		}

		// Store relevant information for interpolation at any timepoint later
		RetrieveCVodeInterpolationInfo();
		cvode_steps++;
		min_step_size = (std::min)(min_step_size, current_simulation_time - prev_time);

		if (DNA_replication_ix != std::numeric_limits<size_t>::max() && replication_start_time != replication_start_time) {
			Real s = NV_Ith_S(cvode_y, DNA_replication_ix);
			if (s > 1e-4) {
				replication_start_time = current_simulation_time;
			}
		}
		if (DNA_replicated_ix != std::numeric_limits<size_t>::max() && replication_finish_time != replication_finish_time) {
			Real s = NV_Ith_S(cvode_y, DNA_replicated_ix);
			if (s > 1.95) {
				replication_finish_time = current_simulation_time;
			}
		}
		if (PCNA_gfp_ix != std::numeric_limits<size_t>::max() && PCNA_gfp_increase_time != PCNA_gfp_increase_time) {
			Real s = NV_Ith_S(cvode_y, PCNA_gfp_ix);
			if (s > 0.5) {
				PCNA_gfp_increase_time = current_simulation_time;
			}
		}
		if (nuclear_envelope_ix != std::numeric_limits<size_t>::max() && nuclear_envelope_breakdown_time != nuclear_envelope_breakdown_time) {
			Real s = NV_Ith_S(cvode_y, nuclear_envelope_ix);
			if (s < 0.5) {
				nuclear_envelope_breakdown_time = current_simulation_time;
			}
		}
		if (chromatid_separation_ix != std::numeric_limits<size_t>::max() && anaphase_onset_time != anaphase_onset_time) {
			if (NV_Ith_S(cvode_y, chromatid_separation_ix) > 1e-3) {
				anaphase_onset_time = current_simulation_time;
			}
		}
		if (cytokinesis_ix != std::numeric_limits<size_t>::max()) {
			if (NV_Ith_S(cvode_y, cytokinesis_ix) > 1.0) {
				divide = true;
			}
		}
		if (apoptosis_ix != std::numeric_limits<size_t>::max()) {
			if (NV_Ith_S(cvode_y, apoptosis_ix) > 1.0) {
				die = true;
			}
		}

		if (divide || die) {
			completed = true;
			break;
		}
	}

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

bool Cell::CellAliveAtTime(Real time, ESynchronizeCellTrajectory synchronize)
{
	if (cvode_steps == 0) {
		return false;
	}
	Real cell_time = time - creation_time;
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

void Cell::SetTreatmentConcentration(Real t)
{
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
		cvode_timepoints_zn[j].col(cvode_steps) = *NV_CONTENT_S(cv_mem->cv_zn[j]);
	}
}
