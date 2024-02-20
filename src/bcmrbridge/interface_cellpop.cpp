#include "Utils.h"
#include "CellPopulationLikelihood.h"
#include "DataLikelihoodDuration.h"
#include "DataLikelihoodTimeCourse.h"
#include "interface.h"

std::shared_ptr<CellPopulationLikelihood> GetCellPopulationLikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<CellPopulationLikelihood>();
	}

	std::shared_ptr<CellPopulationLikelihood> ll = std::dynamic_pointer_cast<CellPopulationLikelihood>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a cell population likelihood");
		*retval = -2;
		return std::shared_ptr<CellPopulationLikelihood>();
	}
	return ll;
}

bool EvaluateCellPopulationLikelihood(bcm3info* info, std::shared_ptr<CellPopulationLikelihood> ll, double* param_values, double* logl, int* retval)
{
	VectorReal param_vector(info->varset->GetNumVariables());
	for (size_t i = 0; i < info->varset->GetNumVariables(); i++) {
		param_vector(i) = param_values[i];
	}

	Real logp = std::numeric_limits<Real>::quiet_NaN();
	if (!ll->EvaluateLogProbability(0, param_vector, logp)) {
		*retval = -3;
		return false;
	}
	if (logl) {
		*logl = logp;
	}

	return true;
}

extern "C" {

void bcm3_rbridge_cellpop_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* logl, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluateCellPopulationLikelihood(info, ll, param_values, logl, retval)) {
		*retval = -4;
		return;
	}
	//ll->OutputEvaluationStatistics("D:");
	*retval = 0;
}

void bcm3_rbridge_cellpop_get_num_species(char** bcm3info_ptr, char** experiment, int* out_num_species, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}
	size_t num_sim_species = e->GetNumSpecies();
	*out_num_species = (int)num_sim_species;
	*retval = 0;
}

void bcm3_rbridge_cellpop_get_species_name(char** bcm3info_ptr, char** experiment, int* species_ix, char** out_name, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}
	snprintf(*out_name, 128, "%s", e->GetSpeciesName(*species_ix).c_str());
	*retval = 0;
}

void bcm3_rbridge_cellpop_get_simulated_trajectories(char** bcm3info_ptr, char** experiment, double* param_values, double* out_values, double* out_timepoints, int* out_parents, int* out_num_cells, int* out_num_timepoints, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}

	Real logl;
	if (!EvaluateCellPopulationLikelihood(info, ll, param_values, &logl, retval)) {
		*retval = -4;
		return;
	}

	size_t num_sim_species = e->GetNumSpecies();

	size_t numcells = e->GetOutputNumCells();
	if (numcells > *out_num_cells) {
		// Need to increase size in the R buffer
		*retval = -5;
		return;
	}
	*out_num_cells = numcells;

	const VectorReal& timepoints = e->GetOutputTimepoints();
	if (timepoints.size() > *out_num_timepoints) {
		// Need to increase size in the R buffer
		*retval = -6;
		return;
	}
	*out_num_timepoints = timepoints.size();
	for (size_t time_i = 0; time_i < timepoints.size(); time_i++) {
		out_timepoints[time_i] = timepoints(time_i);
	}

	for (size_t cell_i = 0; cell_i < numcells; cell_i++) {
		for (size_t time_i = 0; time_i < timepoints.size(); time_i++) {
			const VectorReal& x = e->GetSimulatedTrajectory(cell_i, time_i);
			for (size_t j = 0; j < num_sim_species; j++) {
				out_values[cell_i * timepoints.size() * num_sim_species + time_i * num_sim_species + j] = x(j);
			}
		}
		out_parents[cell_i] = e->GetSimulatedCellParent(cell_i);
	}

	*retval = 0;
}

void bcm3_rbridge_cellpop_get_num_data(char** bcm3info_ptr, char** experiment, int* out_num_data, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}
	size_t ndata = e->GetNumData();
	*out_num_data = ndata;
	*retval = 0;
}

void bcm3_rbridge_cellpop_get_observed_data(char** bcm3info_ptr, char** experiment, int* data_ix, double* out_values, double* out_timepoints, int* out_num_cells, int* out_num_timepoints, int* out_num_markers, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}

	const DataLikelihoodBase* dl = e->GetData(*data_ix);
	const DataLikelihoodTimeCourse* dltc = dynamic_cast<const DataLikelihoodTimeCourse*>(dl);
	if (dltc != nullptr) {
		const VectorReal& t = dltc->GetTimepoints();
		if (t.size() > *out_num_timepoints) {
			// Need to increase size in the R buffer
			*retval = -6;
			return;
		}
		*out_num_timepoints = t.size();
		for (size_t i = 0; i < t.size(); i++) {
			out_timepoints[i] = t(i);
		}

		if (dltc->GetNumObservedData() > *out_num_cells) {
			// Need to increase size in the R buffer
			*retval = -5;
			return;
		}
		*out_num_cells = dltc->GetNumObservedData();

		for (int i = 0; i < *out_num_cells; i++) {
			const MatrixReal& o = dltc->GetObservedData(i);
			*out_num_markers = o.cols();
			for (int j = 0; j < *out_num_markers; j++) {
				for (int k = 0; k < *out_num_timepoints; k++) {
					out_values[i * (*out_num_markers) * (*out_num_timepoints) + j * (*out_num_timepoints) + k] = o(k, j);
				}
			}
		}
	}

	const DataLikelihoodDuration* dld = dynamic_cast<const DataLikelihoodDuration*>(dl);
	if (dld != nullptr) {
		*out_num_timepoints = 1;
		*out_num_cells = dld->GetObservedData().size();
		*out_num_markers = 1;
		for (size_t i = 0; i < *out_num_cells; i++) {
			out_values[i] = dld->GetObservedData()(i);
		}
	}
}

void bcm3_rbridge_cellpop_get_simulated_data(char** bcm3info_ptr, char** experiment, double* param_values, int* data_ix, double* out_values, double* out_timepoints, int* out_num_cells, int* out_num_timepoints, int* out_num_markers, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}

	Real logl;
	if (!EvaluateCellPopulationLikelihood(info, ll, param_values, &logl, retval)) {
		*retval = -4;
		return;
	}

	const DataLikelihoodBase* dl = e->GetData(*data_ix);
	const DataLikelihoodTimeCourse* dltc = dynamic_cast<const DataLikelihoodTimeCourse*>(dl);
	if (dltc != nullptr) {
		const VectorReal& t = dltc->GetTimepoints();
		if (t.size() > *out_num_timepoints) {
			// Need to increase size in the R buffer
			*retval = -6;
			return;
		}
		*out_num_timepoints = t.size();
		for (size_t i = 0; i < t.size(); i++) {
			out_timepoints[i] = t(i);
		}

		if (dltc->GetNumObservedData() > *out_num_cells) {
			// Need to increase size in the R buffer
			*retval = -5;
			return;
		}
		*out_num_cells = dltc->GetNumObservedData();
		for (int i = 0; i < *out_num_cells; i++) {
			const MatrixReal& o = dltc->GetSimulatedData(i);
			*out_num_markers = o.cols();
			for (int j = 0; j < *out_num_markers; j++) {
				for (int k = 0; k < *out_num_timepoints; k++) {
					out_values[i * (*out_num_markers) * (*out_num_timepoints) + j * (*out_num_timepoints) + k] = o(k, j);
				}
			}
		}
	}

	const DataLikelihoodDuration* dld = dynamic_cast<const DataLikelihoodDuration*>(dl);
	if (dld != nullptr) {
		*out_num_timepoints = 1;
		*out_num_cells = dld->GetSimulatedData().size();
		*out_num_markers = 1;
		for (size_t i = 0; i < *out_num_cells; i++) {
			out_values[i] = dld->GetSimulatedData()(i);
		}
	}
}

void bcm3_rbridge_cellpop_get_matched_simulation(char** bcm3info_ptr, char** experiment, double* param_values, int* data_ix, double* out_values, double* out_timepoints, int* out_num_cells, int* out_num_timepoints, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<CellPopulationLikelihood> ll = GetCellPopulationLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const Experiment* e = ll->GetExperiment(*experiment);
	if (!e) {
		*retval = -3;
		return;
	}

	Real logl;
	if (!EvaluateCellPopulationLikelihood(info, ll, param_values, &logl, retval)) {
		*retval = -4;
		return;
	}

	const VectorReal& timepoints = e->GetOutputTimepoints();
	if (timepoints.size() > *out_num_timepoints) {
		// Need to increase size in the R buffer
		*retval = -6;
		return;
	}
	*out_num_timepoints = timepoints.size();
	for (size_t time_i = 0; time_i < timepoints.size(); time_i++) {
		out_timepoints[time_i] = timepoints(time_i);
	}

	const DataLikelihoodBase* dl = e->GetData(*data_ix);
	const DataLikelihoodTimeCourse* dltc = dynamic_cast<const DataLikelihoodTimeCourse*>(dl);
	if (dltc != nullptr) {
		if (dltc->GetNumObservedData() > *out_num_cells) {
			// Need to increase size in the R buffer
			*retval = -5;
			return;
		}
		*out_num_cells = dltc->GetNumObservedData();
		for (int i = 0; i < *out_num_cells; i++) {
			size_t matched_ix = dltc->GetTrajectoryMatching()[i];
			for (size_t time_i = 0; time_i < timepoints.size(); time_i++) {
				const VectorReal& x = e->GetSimulatedTrajectory(matched_ix, time_i);
				for (size_t j = 0; j < e->GetNumSpecies(); j++) {
					out_values[i * timepoints.size() * e->GetNumSpecies() + time_i * e->GetNumSpecies() + j] = x(j);
				}
			}
		}
	}
#if 0
	// TODO
	const DataLikelihoodDuration* dld = dynamic_cast<const DataLikelihoodDuration*>(dl);
	if (dld != nullptr) {
		*out_num_timepoints = 1;
		*out_num_cells = dld->GetSimulatedData().size();
		*out_num_markers = 1;
		for (size_t i = 0; i < *out_num_cells; i++) {
			out_values[i] = dld->GetSimulatedData()(i);
		}
	}
#endif
}

}
