#include "Utils.h"
#include "fISALikelihood.h"
#include "fISAExperiment.h"
#include "interface.h"

std::shared_ptr<fISALikelihood> GetfISALikelihood(bcm3info* info, int* retval)
{
	if (!info) {
		return std::shared_ptr<fISALikelihood>();
	}

	std::shared_ptr<fISALikelihood> ll = std::dynamic_pointer_cast<fISALikelihood>(info->likelihood);
	if (!ll) {
		LOGERROR("bcm3 info pointer does not contain a dynamicISA likelihood");
		*retval = -2;
		return std::shared_ptr<fISALikelihood>();
	}
	return ll;
}

bool EvaluatefISALikelihood(bcm3info* info, std::shared_ptr<fISALikelihood> ll, double* param_values, double* logl, int* retval)
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

void bcm3_rbridge_fISA_get_num_data(char** bcm3info_ptr, char** experiment, int* out_num_data, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const fISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -4;
		return;
	}
	*out_num_data = (int)exp->GetNumData();
	*retval = 0;
}

void bcm3_rbridge_fISA_get_num_cell_lines(char** bcm3info_ptr, char** experiment, int* out_num_cell_lines, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const fISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	*out_num_cell_lines = (int)exp->GetNumCellLines();
	*retval = 0;
}

void bcm3_rbridge_fISA_get_num_replicates(char** bcm3info_ptr, char** experiment, int* data_ix, int* out_num_replicates, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	const fISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	*out_num_replicates = (int)exp->GetNumReplicates(*data_ix);
	*retval = 0;
}

void bcm3_rbridge_fISA_get_num_signaling_molecules(char** bcm3info_ptr, int* out_num_signaling_molecules, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	*out_num_signaling_molecules = (int)ll->GetNumSignalingMolecules();
	*retval = 0;
}

void bcm3_rbridge_fISA_get_cell_line_name(char** bcm3info_ptr, char** experiment, int* cell_line_ix, char** out_cell_line_name, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}
	const fISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	size_t ix = *cell_line_ix;
	if (ix >= exp->GetNumCellLines()) {
		*retval = -4;
		return;
	}
	const std::string& name = exp->GetCellLineName(ix);
	//strncpy_s(*out_cell_line_name, 128, name.c_str(), 128);
	snprintf(*out_cell_line_name, 128, "%s", name.c_str());

	*retval = 0;
}

void bcm3_rbridge_fISA_get_observed_data(char** bcm3info_ptr, char** experiment, int* data_ix, double* out_observed_data_values, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	MatrixReal values;
	const fISAExperiment* exp = ll->GetExperiment(*experiment);
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

void bcm3_rbridge_fISA_get_log_likelihood(char** bcm3info_ptr, double* param_values, double* out_logl, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}
	if (!EvaluatefISALikelihood(info, ll, param_values, out_logl, retval)) {
		*retval = -4;
		return;
	}
	*retval = 0;
}

void bcm3_rbridge_fISA_get_modeled_data(char** bcm3info_ptr, char** experiment, double* param_values, int* data_ix, double* out_modeled_data_values, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluatefISALikelihood(info, ll, param_values, nullptr, retval)) {
		*retval = -4;
		return;
	}

	VectorReal values;
	const fISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	exp->GetModeledData(0, (size_t)*data_ix, values);
	for (size_t i = 0; i < values.size(); i++) {
		out_modeled_data_values[i] = values(i);
	}

	*retval = 0;
}

void bcm3_rbridge_fISA_get_modeled_activities(char** bcm3info_ptr, char** experiment, double* param_values, double* out_modeled_activities, int* retval)
{
	bcm3info* info = GetBCM3InfoPtr(bcm3info_ptr, retval);
	std::shared_ptr<fISALikelihood> ll = GetfISALikelihood(info, retval);
	if (!info || !ll) {
		return;
	}

	if (!EvaluatefISALikelihood(info, ll, param_values, nullptr, retval)) {
		*retval = -4;
		return;
	}

	MatrixReal values;
	const fISAExperiment* exp = ll->GetExperiment(*experiment);
	if (!exp) {
		*retval = -3;
		return;
	}
	exp->GetModeledActivities(0, values);
	for (size_t i = 0; i < values.rows(); i++) {
		for (size_t j = 0; j < values.cols(); j++) {
			out_modeled_activities[i * values.cols() + j] = values(i, j);
		}
	}

	*retval = 0;
}

}