#include "Utils.h"
#include "Experiment.h"
#include "TreatmentTrajectoryPulses.h"

bool TreatmentTrajectoryPulses::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::NetCDFDataFile& data_file)
{
	std::string times_str = xml_node.get<std::string>("<xmlattr>.times");

	std::vector<std::string> time_tokens;
	bcm3::tokenize(times_str, time_tokens, ",");
	timepoints.setConstant(time_tokens.size(), std::numeric_limits<Real>::quiet_NaN());
	for (int i = 0; i < time_tokens.size(); i++) {
		timepoints(i) = boost::lexical_cast<Real>(time_tokens[i]);
	}

	std::sort(timepoints.begin(), timepoints.end());

	return true;
}

Real TreatmentTrajectoryPulses::GetConcentration(Real time, Real cell_creation_time)
{
	Real global_time = time + cell_creation_time;

	for (int i = 0; i < timepoints.size(); i++) {
		Real t_in_pulse = global_time - timepoints(i) - 2.0;
		if (t_in_pulse >= 14.0 || t_in_pulse <= 0.0) {
			return 0.0;
		} else if (t_in_pulse < 2.0) {
			return t_in_pulse * 0.5;
		} else if (t_in_pulse < 10.0) {
			return 1.0;
		} else {
			return 1 - (t_in_pulse - 10.0) * 0.25;
		}
	}

	return 0.0;
}

Real TreatmentTrajectoryPulses::FirstDiscontinuity(Real cell_creation_time)
{
	if (timepoints.size() > 0) {
		return timepoints[0] - cell_creation_time + 2.0;
	} else {
		return std::numeric_limits<Real>::quiet_NaN();
	}
}

Real TreatmentTrajectoryPulses::NextDiscontinuity(Real time, Real cell_creation_time)
{
	for (int i = 0; i < timepoints.size(); i++) {
		if (time == timepoints[i] - cell_creation_time + 2.0) {
			return timepoints[i] - cell_creation_time + 4.0;
		} else if (time == timepoints[i] - cell_creation_time + 4.0) {
			return timepoints[i] - cell_creation_time + 10.0;
		} else if (time == timepoints[i] - cell_creation_time + 10.0) {
			return timepoints[i] - cell_creation_time + 14.0;
		} else if (time == timepoints[i] - cell_creation_time + 14.0) {
			if (i < timepoints.size() - 1) {
				return timepoints[i + 1] - cell_creation_time + 2.0;
			} else {
				return std::numeric_limits<OdeReal>::quiet_NaN();
			}
		}
	}

	return std::numeric_limits<OdeReal>::quiet_NaN();
}
