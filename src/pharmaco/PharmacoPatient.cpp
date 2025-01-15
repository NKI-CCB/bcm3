#include "Utils.h"
#include "PharmacoPatient.h"

Patient::Patient()
{
}

bool Patient::Load(const bcm3::NetCDFDataFile& data, const std::string& trial, const std::string& drug)
{
	bool result = true;

	// Find patient index
	size_t patient_ix;
	if (!data.GetDimensionIx(trial, "patients", patient_id, &patient_ix)) {
		LOGERROR("Cannot find patient \"%s\" in data file", patient_id.c_str());
		return false;
	}

	Real dose;
	Real dosing_interval;
	Real dose_after_dose_change;
	Real dose_change_time;
	unsigned int intermittent;
	std::set<int> skipped_days;

	size_t num_timepoints;
	result &= data.GetDimensionSize(trial, "time", &num_timepoints);

	observation_timepoints.resize(num_timepoints);
	observed_concentrations.resize(num_timepoints);
	result &= data.GetValues(trial, "time", 0, num_timepoints, observation_timepoints);
	result &= data.GetValuesDim2(trial, drug + "_plasma_concentration", patient_ix, 0, num_timepoints, observed_concentrations);
	result &= data.GetValue(trial, drug + "_dose", patient_ix, &dose);
	result &= data.GetValue(trial, drug + "_dosing_interval", patient_ix, &dosing_interval);
	result &= data.GetValue(trial, drug + "_dose_after_dose_change", patient_ix, &dose_after_dose_change);
	result &= data.GetValue(trial, drug + "_dose_change_time", patient_ix, &dose_change_time);
	result &= data.GetValue(trial, drug + "_intermittent", patient_ix, &intermittent);

	std::vector<unsigned int> interruptions;
	for (int i = 0; i < 29; i++) {
		unsigned int tmp;
		result &= data.GetValue(trial, "treatment_interruptions", patient_ix, (size_t)i, &tmp);
		if (tmp) {
			skipped_days.insert(i);
		}
	}

	std::vector<Real> times;
	Real last_time = 696;
	Real t = 0;
	while (t < last_time) {
		bool give_treatment = true;
		int day = static_cast<int>(floor(t / 24.0));
		if (skipped_days.find(day) != skipped_days.end()) {
			give_treatment = false;
		}
		if (intermittent == 1) {
			Real time_in_week = t - 7.0 * 24.0 * floor(t / (7.0 * 24.0));
			if (time_in_week >= 5.0 * 24.0) {
				// No dose in day 6 and 7
				give_treatment = false;
			}
		} else if (intermittent == 2) {
			Real time_in_treatment_course = t - 28.0 * 24.0 * floor(t / (28.0 * 24.0));
			if (time_in_treatment_course >= 21.0 * 24.0) {
				// No dose in day 22 through 28
				give_treatment = false;
			}
		} else if (intermittent == 3) {
			Real time_in_week = t - 7.0 * 24.0 * floor(t / (7.0 * 24.0));
			if (time_in_week >= 4.0 * 24.0) {
				// No dose in day 5, 6 and 7
				give_treatment = false;
			}
		}

		if (give_treatment) {
			times.push_back(t);
		}
		t += dosing_interval;
	}
	treatment_timepoints.resize(times.size());
	treatment_doses.resize(times.size());
	for (ptrdiff_t i = 0; i < times.size(); i++) {
		treatment_timepoints(i) = times[i];
		if (!std::isnan(dose_change_time) && times[i] >= dose_change_time) {
			treatment_doses(i) = dose_after_dose_change;
		} else {
			treatment_doses(i) = dose;
		}
	}

	VectorReal new_observed_timepoints = observation_timepoints;
	VectorReal new_observed_concentrations = observed_concentrations;
	ptrdiff_t out_i = 0;
	for (ptrdiff_t i = 0; i < observation_timepoints.size(); i++) {
		if (!std::isnan(observed_concentrations(i))) {
			new_observed_timepoints(out_i) = observation_timepoints(i);
			new_observed_concentrations(out_i) = observed_concentrations(i);
			out_i++;
		}
	}
	observation_timepoints = new_observed_timepoints.segment(0, out_i);
	observed_concentrations = new_observed_concentrations.segment(0, out_i);
	simulated_concentrations.resize(observed_concentrations.size());

	return true;
}
