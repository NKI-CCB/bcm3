#include "Utils.h"
#include "LikelihoodODE.h"
#include "interface.h"

std::shared_ptr<LikelihoodODE> GetODELikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<LikelihoodODE>();
	}

	std::shared_ptr<LikelihoodODE> ll = std::dynamic_pointer_cast<LikelihoodODE>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a cell population likelihood");
		*retval = -2;
		return std::shared_ptr<LikelihoodODE>();
	}
	return ll;
}

bool EvaluateLikelihoodODE(bcm3info* info, std::shared_ptr<LikelihoodODE> ll, double* param_values, double* logl, int* retval)
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

	void bcm3_rbridge_ODE_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* logl, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<LikelihoodODE> ll = GetODELikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		if (!EvaluateLikelihoodODE(info, ll, param_values, logl, retval)) {
			*retval = -4;
			return;
		}
		*retval = 0;
	}

	void bcm3_rbridge_ODE_get_simulated_trajectories(char** bcm3info_ptr, double* param_values, double* out_values, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<LikelihoodODE> ll = GetODELikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		Real logl;
		if (!EvaluateLikelihoodODE(info, ll, param_values, &logl, retval)) {
			*retval = -4;
			return;
		}

		const OdeMatrixReal& simtraj = ll->GetSimulatedTrajectories();
		for (size_t i = 0; i < 100; i++) {
			for (size_t j = 0; j < 4; j++) {
				out_values[j * 100 + i] = simtraj(j, i);
			}
		}

		*retval = 0;
	}

}
