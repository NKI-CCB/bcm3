#include "Utils.h"
#include "PharmacoLikelihoodSingle.h"
#include "interface.h"

std::shared_ptr<PharmacoLikelihoodSingle> GetPKLikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<PharmacoLikelihoodSingle>();
	}

	std::shared_ptr<PharmacoLikelihoodSingle> ll = std::dynamic_pointer_cast<PharmacoLikelihoodSingle>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a population PK likelihood");
		*retval = -2;
		return std::shared_ptr<PharmacoLikelihoodSingle>();
	}
	return ll;
}

bool EvaluatePKLikelihood(bcm3info* info, std::shared_ptr<PharmacoLikelihoodSingle> ll, double* param_values, double* logl, int* retval)
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

	void bcm3_rbridge_pharmacosingle_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* logl, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<PharmacoLikelihoodSingle> ll = GetPKLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		if (!EvaluatePKLikelihood(info, ll, param_values, logl, retval)) {
			*retval = -4;
			return;
		}
		*retval = 0;
	}

	void bcm3_rbridge_pharmacosingle_get_observed_data(char** bcm3info_ptr, double* out_values, double* out_timepoints, int* out_num_timepoints, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<PharmacoLikelihoodSingle> ll = GetPKLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		const VectorReal& t = ll->GetTimepoints();
		if (t.size() > *out_num_timepoints) {
			LOGERROR("Insufficient space in output buffer");
			*retval = -5;
			return;
		}
		*out_num_timepoints = t.size();
		for (ptrdiff_t i = 0; i < t.size(); i++) {
			out_timepoints[i] = t(i);
		}

		for (ptrdiff_t i = 0; i < t.size(); i++) {
			out_values[i] = ll->GetObservedConcentrations()(i);
		}
		*retval = 0;
	}

	void bcm3_rbridge_pharmacosingle_get_simulated_data(char** bcm3info_ptr, double* param_values, double* out_values, double* out_timepoints, int* out_num_timepoints, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<PharmacoLikelihoodSingle> ll = GetPKLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}
		Real logl;
		if (!EvaluatePKLikelihood(info, ll, param_values, &logl, retval)) {
			*retval = -4;
			return;
		}

		const VectorReal& t = ll->GetTimepoints();
		if (t.size() > *out_num_timepoints) {
			LOGERROR("Insufficient space in output buffer");
			*retval = -5;
			return;
		}
		*out_num_timepoints = t.size();
		for (ptrdiff_t i = 0; i < t.size(); i++) {
			out_timepoints[i] = t(i);
		}

		for (ptrdiff_t i = 0; i < t.size(); i++) {
			out_values[i] = ll->GetSimulatedConcentrations()(i);
		}
		*retval = 0;
	}

	void bcm3_rbridge_pharmacosingle_get_simulated_trajectory(char** bcm3info_ptr, double* param_values, double* timepoints, int* num_timepoints, double* out_trajectories, double* out_concentrations, int* out_num_compartments, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<PharmacoLikelihoodSingle> ll = GetPKLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}
		Real logl;
		if (!EvaluatePKLikelihood(info, ll, param_values, &logl, retval)) {
			*retval = -4;
			return;
		}

		VectorReal time(*num_timepoints);
		for (int i = 0; i < *num_timepoints; i++) {
			time(i) = timepoints[i];
		}

		VectorReal concentrations(*num_timepoints);
		MatrixReal trajectory;
		if (!ll->GetSimulatedTrajectory(time, concentrations, trajectory)) {
			*retval = -6;
			return;
		}

		if (trajectory.cols() > *out_num_compartments) {
			LOGERROR("Insufficient compartments in output buffer");
			*retval = -5;
			return;
		}
		*out_num_compartments = trajectory.cols();

		for (int i = 0; i < *num_timepoints; i++) {
			const Real t = timepoints[i];
			out_concentrations[i] = concentrations(i);

			for (int j = 0; j < *out_num_compartments; j++) {
				out_trajectories[i * (*out_num_compartments) + j] = trajectory(i, j);
			}
		}

		*retval = 0;
	}

}
