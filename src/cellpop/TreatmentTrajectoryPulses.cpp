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

	// Iterate backwards because the pulse calculation returns most quickly if time is before the pulse time
	Real value = 0.0;
	for (int i = timepoints.size() - 1; i >= 0; i--) {
		value += CalculatePulse(global_time, timepoints(i));
		if (value > 0.0) {
			// Pulses do not overlap so we can quickly bail out
			break;
		}
	}

	return value;
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

Real TreatmentTrajectoryPulses::CalculatePulse(Real t, Real pulse_time)
{
	Real t_in_pulse = t - pulse_time - 2.0;

	//return ((0.5 * tanh(100.0 * t_in - 100.0 * t + 300.0) + 0.5) * (t_in - 1.0 * t + 2.0) + (0.5 * tanh(100.0 * t_in - 100.0 * t + 300.0) - 0.5) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1000.0) - 1.0 * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1400.0) + 0.5) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 1000.0) - 0.5) * (0.25 * t_in - 0.25 * t + 3.5) + 0.5)) * (0.5 * tanh(100.0 * t_in - 100.0 * t + 200.0) - 0.5);
	if (t_in_pulse < 0.0) {
		return 0.0;
	} else if (t_in_pulse < 2.0) {
		return t_in_pulse * 0.5;
	} else if (t_in_pulse < 10.0) {
		return 1.0;
	} else if (t_in_pulse < 14.0) {
		return 1 - (t_in_pulse - 10.0) * 0.25;
	} else {
		return 0.0;
	}
}
