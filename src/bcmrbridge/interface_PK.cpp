#include "Utils.h"
#include "LikelihoodPharmacokineticTrajectory.h"
#include "interface.h"

std::shared_ptr<LikelihoodPharmacokineticTrajectory> GetPopPKLikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<LikelihoodPharmacokineticTrajectory>();
	}

	std::shared_ptr<LikelihoodPharmacokineticTrajectory> ll = std::dynamic_pointer_cast<LikelihoodPharmacokineticTrajectory>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a population PK likelihood");
		*retval = -2;
		return std::shared_ptr<LikelihoodPharmacokineticTrajectory>();
	}
	return ll;
}

bool EvaluatePKLikelihood(bcm3info* info, std::shared_ptr<LikelihoodPharmacokineticTrajectory> ll, double* param_values, double* logl, int* retval)
{
	VectorReal param_vector(info->varset->GetNumVariables());
	for (size_t i = 0; i < info->varset->GetNumVariables(); i++) {
		param_vector(i) = param_values[i];
	}

	Real logp;
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

void bcm3_rbridge_PK_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* logl, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<LikelihoodPharmacokineticTrajectory> ll = GetPopPKLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluatePKLikelihood(info, ll, param_values, logl, retval)) {
		*retval = -4;
		return;
	}
	*retval = 0;
}

void bcm3_rbridge_PK_get_observed_data(char** bcm3info_ptr, double* out_values, double* out_timepoints, int* out_num_timepoints, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<LikelihoodPharmacokineticTrajectory> ll = GetPopPKLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const OdeVectorReal& t = ll->GetTimepoints();
	*out_num_timepoints = t.size();
	for (ptrdiff_t i = 0; i < t.size(); i++) {
		out_timepoints[i] = t(i);
	}

	for (ptrdiff_t i = 0; i < t.size(); i++) {
		out_values[i] = ll->GetObservedConcentrations()(i);
	}
	*retval = 0;
}

void bcm3_rbridge_PK_get_simulated_data(char** bcm3info_ptr, double* param_values, double* out_trajectories, double* out_concentrations, double* out_timepoints, int* out_num_timepoints, int* out_num_compartments, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<LikelihoodPharmacokineticTrajectory> ll = GetPopPKLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}
	Real logl;
	if (!EvaluatePKLikelihood(info, ll, param_values, &logl, retval)) {
		*retval = -4;
		return;
	}

	const OdeVectorReal& t = ll->GetTimepoints();
	*out_num_timepoints = t.size();
	for (ptrdiff_t i = 0; i < t.size(); i++) {
		out_timepoints[i] = t(i);
	}

	*out_num_compartments = (int)ll->GetSimulatedTrajectories(0).rows();
	const OdeVectorReal& simulated_concentrations = ll->GetSimulatedConcentrations(0);
	const OdeMatrixReal& simulated_trajectories = ll->GetSimulatedTrajectories(0);
	for (ptrdiff_t i = 0; i < simulated_concentrations.size(); i++) {
		out_concentrations[i] = simulated_concentrations(i);
		for (int k = 0; k < (*out_num_compartments); k++) {
			out_trajectories[k * t.size() + i] = simulated_trajectories(k, i);
		}
	}
	if (simulated_concentrations.size() < t.size()) {
		for (ptrdiff_t i = simulated_concentrations.size(); i < t.size(); i++) {
			out_concentrations[i] = std::numeric_limits<double>::quiet_NaN();
			for (int k = 0; k < (*out_num_compartments); k++) {
				out_trajectories[k * t.size() + i] = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}

	*retval = 0;
}

}