#include "Utils.h"
#include "Experiment.h"
#include "TreatmentTrajectory.h"

bool TreatmentTrajectory::Load(Experiment* experiment, const bcm3::NetCDFDataFile& data_file, const std::string& treatment_variable)
{
	bool result = true;

	size_t num_timepoints, num_samples;
	result &= data_file.GetDimensionSize(experiment->GetName(), "treatment_time", &num_timepoints);
	result &= data_file.GetDimensionSize(experiment->GetName(), "treatment_samples", &num_samples);

	timepoints.resize(num_timepoints);
	concentrations.resize(num_samples, num_timepoints);

	data_file.GetValues(experiment->GetName(), "treatment_time", 0, num_timepoints, timepoints);
	timepoints *= 3600;

	for (size_t i = 0; i < num_samples; i++) {
		VectorReal v;
		data_file.GetValuesDim2(experiment->GetName(), treatment_variable, i, 0, num_timepoints, v);
		concentrations.row(i) = v;
	}

	return result;
}

Real TreatmentTrajectory::GetConcentration(Real time, size_t sample_ix)
{
	ASSERT(timepoints.size() >= 2);
	ASSERT(timepoints.size() == concentrations.cols());
	ASSERT(sample_ix < concentrations.rows());

	if (time <= timepoints(0)) {
		return concentrations(sample_ix, 0);
	} else if (time >= timepoints(timepoints.size() - 1)) {
		return concentrations(sample_ix, timepoints.size() - 1);
	}

	for (size_t i = 0; i < timepoints.size()-1; i++) {
		if (timepoints(i) == time) {
			return concentrations(sample_ix, i);
		} else if (timepoints(i) < time && timepoints(i+1) > time) {
			Real x = (time - timepoints(i)) / (timepoints(i+1) - timepoints(i));
			return concentrations(sample_ix, i) + x * (concentrations(sample_ix, i+1) - concentrations(sample_ix, i));
		}
	}

	return std::numeric_limits<Real>::quiet_NaN();
}
