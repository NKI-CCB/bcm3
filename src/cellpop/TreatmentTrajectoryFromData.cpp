#include "Utils.h"
#include "Experiment.h"
#include "TreatmentTrajectoryFromData.h"

bool TreatmentTrajectoryFromData::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::NetCDFDataFile& data_file)
{
	if (!TreatmentTrajectory::Load(xml_node, experiment, data_file)) {
		return false;
	}

	bool result = true;

	size_t num_timepoints, num_samples;
	result &= data_file.GetDimensionSize(experiment->GetName(), "treatment_time", &num_timepoints);
	result &= data_file.GetDimensionSize(experiment->GetName(), "treatment_samples", &num_samples);

	timepoints.resize(num_timepoints);
	concentrations.resize(num_samples, num_timepoints);

	data_file.GetValues(experiment->GetName(), "treatment_time", 0, num_timepoints, timepoints);
	timepoints *= 3600;

#if TODO
	// Load treatment_variable from xml_node

	for (size_t i = 0; i < num_samples; i++) {
		VectorReal v;
		data_file.GetValuesDim2(experiment->GetName(), treatment_variable, i, 0, num_timepoints, v);
		concentrations.row(i) = v;
	}
	return result;
#else
	return false;
#endif
}

Real TreatmentTrajectoryFromData::GetConcentration(Real time, Real cell_creation_time)
{
	Real global_time = time + cell_creation_time;

	ASSERT(timepoints.size() >= 2);
	ASSERT(timepoints.size() == concentrations.cols());

#if TODO
#else
	size_t sample_ix = 0;
#endif

	if (global_time <= timepoints(0)) {
		return concentrations(sample_ix, 0);
	} else if (global_time >= timepoints(timepoints.size() - 1)) {
		return concentrations(sample_ix, timepoints.size() - 1);
	}

	for (size_t i = 0; i < timepoints.size()-1; i++) {
		if (timepoints(i) == global_time) {
			return concentrations(sample_ix, i);
		} else if (timepoints(i) < global_time && timepoints(i+1) > global_time) {
			Real x = (global_time - timepoints(i)) / (timepoints(i+1) - timepoints(i));
			return concentrations(sample_ix, i) + x * (concentrations(sample_ix, i+1) - concentrations(sample_ix, i));
		}
	}

	return std::numeric_limits<Real>::quiet_NaN();
}
