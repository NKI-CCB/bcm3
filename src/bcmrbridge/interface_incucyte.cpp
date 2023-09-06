#include "Utils.h"
#include "LikelihoodIncucytePopulation.h"
#include "interface.h"

std::shared_ptr<LikelihoodIncucytePopulation> GetLikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<LikelihoodIncucytePopulation>();
	}

	std::shared_ptr<LikelihoodIncucytePopulation> ll = std::dynamic_pointer_cast<LikelihoodIncucytePopulation>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a cell population likelihood");
		*retval = -2;
		return std::shared_ptr<LikelihoodIncucytePopulation>();
	}
	return ll;
}

bool EvaluateLikelihoodIncucyte(bcm3info* info, std::shared_ptr<LikelihoodIncucytePopulation> ll, double* param_values, double* logl, int* retval)
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

	void bcm3_rbridge_incucyte_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* logl, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<LikelihoodIncucytePopulation> ll = GetLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		if (!EvaluateLikelihoodIncucyte(info, ll, param_values, logl, retval)) {
			*retval = -4;
			return;
		}
		*retval = 0;
	}

	void bcm3_rbridge_incucyte_get_simulated_trajectories(char** bcm3info_ptr, double* param_values, double* out_values, int* experiment_ix, int* num_timepoints, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<LikelihoodIncucytePopulation> ll = GetLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		Real logl;
		if (!EvaluateLikelihoodIncucyte(info, ll, param_values, &logl, retval)) {
			*retval = -4;
			return;
		}

		const MatrixReal matrices[5] = {
			ll->GetSimulatedCellCount(*experiment_ix - 1),
			ll->GetSimulatedApoptoticCellCount(*experiment_ix - 1),
			ll->GetSimulatedDebris(*experiment_ix - 1),
			ll->GetSimulatedConfluence(*experiment_ix - 1),
			ll->GetSimulatedApoptosisMarker(*experiment_ix - 1)
		};
		if (matrices[0].rows() * matrices[0].cols() > 100 * 11) {
			LOGERROR("Simulated trajectories too large - increase buffer size");
			*retval = -5;
			return;
		}

		for (ptrdiff_t k = 0; k < 5; k++) {
			for (ptrdiff_t i = 0; i < matrices[0].rows(); i++) {
				for (ptrdiff_t j = 0; j < matrices[0].cols(); j++) {
					out_values[k * matrices[0].rows() * matrices[0].cols() + j * matrices[0].rows() + i] = matrices[k](i, j);
				}
			}
		}
		*num_timepoints = matrices[0].rows();

		*retval = 0;
	}

	void bcm3_rbridge_incucyte_get_simulated_ctb(char** bcm3info_ptr, double* param_values, double* out_values, int* experiment_ix, int* retval)
	{
		bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
		std::shared_ptr<LikelihoodIncucytePopulation> ll = GetLikelihood(info, retval);
		if (!info || !ll) {
			return;
		}

		Real logl;
		if (!EvaluateLikelihoodIncucyte(info, ll, param_values, &logl, retval)) {
			*retval = -4;
			return;
		}

		const VectorReal ctb = ll->GetSimulatedCTB(*experiment_ix-1);
		if (ctb.size() > 9) {
			LOGERROR("Simulated CTB vector too large - increase buffer size");
			*retval = -5;
			return;
		}

		for (ptrdiff_t i = 0; i < ctb.size(); i++) {
			out_values[i] = ctb(i);
		}
		*retval = 0;
	}

}
