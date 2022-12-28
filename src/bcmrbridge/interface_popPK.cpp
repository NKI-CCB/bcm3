#include "Utils.h"
#include "LikelihoodPopPKTrajectory.h"
#include "interface.h"

std::shared_ptr<LikelihoodPopPKTrajectory> GetPopPKLikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<LikelihoodPopPKTrajectory>();
	}

	std::shared_ptr<LikelihoodPopPKTrajectory> ll = std::dynamic_pointer_cast<LikelihoodPopPKTrajectory>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a population PK likelihood");
		*retval = -2;
		return std::shared_ptr<LikelihoodPopPKTrajectory>();
	}
	return ll;
}

bool EvaluatePopPKLikelihood(bcm3info* info, std::shared_ptr<LikelihoodPopPKTrajectory> ll, double* param_values, double* logl, int* retval)
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

void bcm3_rbridge_popPK_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* logl, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<LikelihoodPopPKTrajectory> ll = GetPopPKLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluatePopPKLikelihood(info, ll, param_values, logl, retval)) {
		*retval = -4;
		return;
	}
	*retval = 0;
}

void bcm3_rbridge_popPK_get_observed_data(char** bcm3info_ptr, double* out_values, double* out_timepoints, int* out_num_patients, int* out_num_timepoints, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<LikelihoodPopPKTrajectory> ll = GetPopPKLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const OdeVectorReal& t = ll->GetTimepoints();
	*out_num_timepoints = t.size();
	for (ptrdiff_t i = 0; i < t.size(); i++) {
		out_timepoints[i] = t(i);
	}

	*out_num_patients = (int)ll->GetNumPatients();
	for (size_t j = 0; j < ll->GetNumPatients(); j++) {
		for (ptrdiff_t i = 0; i < t.size(); i++) {
			out_values[j * t.size() + i] = ll->GetObservedConcentrations(j)(i);
		}
	}
	*retval = 0;
}

void bcm3_rbridge_popPK_get_simulated_data(char** bcm3info_ptr, double* param_values, double* out_trajectories, double* out_concentrations, double* out_timepoints, int* out_num_patients, int* out_num_timepoints, int* out_num_compartments, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<LikelihoodPopPKTrajectory> ll = GetPopPKLikelihood(info, retval);
	if (!info || !ll) {
		return;
	}
	Real logl;
	if (!EvaluatePopPKLikelihood(info, ll, param_values, &logl, retval)) {
		*retval = -4;
		return;
	}

	const OdeVectorReal& t = ll->GetTimepoints();
	*out_num_timepoints = t.size();
	for (ptrdiff_t i = 0; i < t.size(); i++) {
		out_timepoints[i] = t(i);
	}

	*out_num_patients = (int)ll->GetNumPatients();
	*out_num_compartments = (int)ll->GetSimulatedTrajectories(0, 0).rows();
	for (size_t j = 0; j < ll->GetNumPatients(); j++) {
		for (ptrdiff_t i = 0; i < t.size(); i++) {
			out_concentrations[j * t.size() + i] = ll->GetSimulatedConcentrations(0, j)(i);
			for (int k = 0; k < (*out_num_compartments); k++) {
				out_trajectories[j * t.size() * (*out_num_compartments) + k * t.size() + i] = ll->GetSimulatedTrajectories(0, j)(k, i);
			}
		}
	}

	*retval = 0;
}

}