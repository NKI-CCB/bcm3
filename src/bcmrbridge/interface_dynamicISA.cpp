#include "Utils.h"
#include "dynamicISALikelihood.h"
#include "dynamicISAExperiment.h"
#include "interface.h"

std::shared_ptr<dynamicISALikelihood> GetDynamicISALikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<dynamicISALikelihood>();
	}

	std::shared_ptr<dynamicISALikelihood> ll = std::dynamic_pointer_cast<dynamicISALikelihood>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a dynamicISA likelihood");
		*retval = -2;
		return std::shared_ptr<dynamicISALikelihood>();
	}
	return ll;
}

bool EvaluateDynamicISALikelihood(bcm3info* info, std::shared_ptr<dynamicISALikelihood> ll, double* param_values, int* retval)
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
	return true;
}

extern "C" {

void bcm3_rbridge_dynamicISA_get_num_data(char** bcm3info_ptr, char** experiment, int* out_num_data, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const dynamicISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -4;
		return;
	}
	*out_num_data = (int)exp->GetNumData();
	*retval = 0;
}

void bcm3_rbridge_dynamicISA_get_num_conditions(char** bcm3info_ptr, char** experiment, int* out_num_conditions, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const dynamicISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	*out_num_conditions = (int)exp->GetNumConditions();
	*retval = 0;
}

void bcm3_rbridge_dynamicISA_get_num_replicates(char** bcm3info_ptr, char** experiment, int* data_ix, int* out_num_replicates, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const dynamicISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	*out_num_replicates = (int)exp->GetNumReplicates(*data_ix);
	*retval = 0;
}

void bcm3_rbridge_dynamicISA_get_num_signaling_molecules(char** bcm3info_ptr, int* out_num_signaling_molecules, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	*out_num_signaling_molecules = (int)ll->GetNumSignalingMolecules();
	*retval = 0;
}

void bcm3_rbridge_dynamicISA_get_observed_data(char** bcm3info_ptr, char** experiment, int* data_ix, double* out_observed_data_values, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	MatrixReal values;
	const dynamicISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	exp->GetObservedData((size_t)*data_ix, values);
	for (size_t i = 0; i < values.rows(); i++) {
		for (size_t j = 0; j < values.cols(); j++) {
			out_observed_data_values[i * values.cols() + j] = values(i, j);
		}
	}

	*retval = 0;
}

void bcm3_rbridge_dynamicISA_get_modeled_data(char** bcm3info_ptr, char** experiment, double* param_values, int* data_ix, double* out_modeled_data_values, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluateDynamicISALikelihood(info, ll, param_values, retval)) {
		*retval = -4;
		return;
	}

	VectorReal values;
	const dynamicISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	exp->GetModeledData((size_t)*data_ix, values);
	for (size_t i = 0; i < values.size(); i++) {
		out_modeled_data_values[i] = values(i);
	}

	*retval = 0;
}

void bcm3_rbridge_dynamicISA_get_modeled_activities(char** bcm3info_ptr, char** experiment, double* param_values, double* out_modeled_activities, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<dynamicISALikelihood> ll = GetDynamicISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluateDynamicISALikelihood(info, ll, param_values, retval)) {
		*retval = -4;
		return;
	}

	MatrixReal values;
	const dynamicISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	exp->GetModeledActivities(values);
	for (size_t i = 0; i < values.rows(); i++) {
		for (size_t j = 0; j < values.cols(); j++) {
			out_modeled_activities[i * values.cols() + j] = values(i, j);
		}
	}

	*retval = 0;
}

}