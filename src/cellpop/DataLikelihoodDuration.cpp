#include "DataLikelihoodDuration.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
#include "../../dependencies/HungarianAlgorithm-master/matching.h"

DataLikelihoodDuration::DataLikelihoodDuration(size_t parallel_evaluations)
{

}

DataLikelihoodDuration::~DataLikelihoodDuration()
{

}

bool DataLikelihoodDuration::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	if (!DataLikelihoodBase::Load(xml_node, experiment, varset, data_file, vm)) {
		return false;
	}

	Real simulation_time = xml_node.get<Real>("<xmlattr>.simulation_time");
	experiment->AddSimulationTimepoints(this, 0.0, 0, std::numeric_limits<size_t>::max(), ESynchronizeCellTrajectory::None);
	experiment->AddSimulationTimepoints(this, simulation_time, 1, std::numeric_limits<size_t>::max(), ESynchronizeCellTrajectory::None);

	std::string periodstr = xml_node.get<std::string>("<xmlattr>.period");
	EPhaseDuration period = EPhaseDuration::None;
	if (periodstr == "G1phase") {
		period = EPhaseDuration::G1phase;
	} else if (periodstr == "Sphase") {
		period = EPhaseDuration::Sphase;
	} else if (periodstr == "G2phase") {
		period = EPhaseDuration::G2phase;
	} else if (periodstr == "NEBD_to_AnaphaseOnset") {
		period = EPhaseDuration::NEBD_to_AnaphaseOnset;
	} else {
		LOGERROR("Unknown duration period \"%s\"", periodstr.c_str());
		return false;
	}

	bool result = true;

	std::string cell_dimension;
	result &= data_file.GetDimensionName(experiment->GetName(), data_name, 0, cell_dimension);
	size_t num_cells;
	result &= data_file.GetDimensionSize(experiment->GetName(), cell_dimension, &num_cells);
	if (!result) {
		return false;
	}

	result &= data_file.GetValues(experiment->GetName(), data_name, 0, num_cells, observed_durations);
	matched_durations.setConstant(observed_durations.size(), std::numeric_limits<Real>::quiet_NaN());

	simulated_durations.setConstant(experiment->GetMaxNumberOfCells(), std::numeric_limits<Real>::quiet_NaN());
	for (size_t i = 0; i < num_cells; i++) {
		experiment->AddRequestedDuration(this, period);
	}
	used_durations.resize(observed_durations.size(), 0);
	if (experiment->GetMaxNumberOfCells() < num_cells) {
		LOGERROR("Simulating less cells than needed for the likelihood in experiment \"%s\"", experiment->GetName().c_str());
		return false;
	}

	return result;
}

bool DataLikelihoodDuration::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	if (!DataLikelihoodBase::PostInitialize(varset, non_sampled_parameter_names)) {
		return false;
	}

	return true;
}

void DataLikelihoodDuration::Reset()
{
	simulated_durations.setConstant(std::numeric_limits<Real>::quiet_NaN());
	matched_durations.setConstant(std::numeric_limits<Real>::quiet_NaN());
	for (int i = 0; i < observed_durations.size(); i++) {
		used_durations[i] = i;
	}
}

bool DataLikelihoodDuration::Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	Real stdev;
	if (stdev_ix != std::numeric_limits<size_t>::max()) {
		stdev = transformed_values[stdev_ix];
	} else if (non_sampled_stdev_ix != std::numeric_limits<size_t>::max()) {
		stdev = non_sampled_parameters[non_sampled_stdev_ix];
	} else {
		stdev = fixed_stdev_value;
	}
	stdev *= stdev_multiplication_factor;
	stdev += 1e-4;

	int n = std::max(observed_durations.size(), observed_durations.size());
	MatrixReal likelihoods(n, n);
	likelihoods.setZero();
	for (int i = 0; i < observed_durations.size(); i++) {
		Real y = observed_durations(i);
		if (y == y) {
			int count = 0;
			for (int j = 0; j < simulated_durations.size(); j++) {
				Real x = simulated_durations(j);
				if (x == x) {
					likelihoods(i, count) = bcm3::LogPdfNormal(y, x, stdev);
					used_durations[count] = j;
					count++;
					if (count == observed_durations.size()) {
						break;
					}
				}
			}
			if (count != observed_durations.size()) {
				logp = -std::numeric_limits<Real>::infinity();
				return true;
			}
		}
	}

	MatrixReal matching;
	if (!HungarianMatching(likelihoods, matching, HUNGARIAN_MATCH_MAX)) {
		LOGERROR("Could not find matching");
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	} else {
		for (int i = 0; i < observed_durations.size(); i++) {
			for (int j = 0; j < observed_durations.size(); j++) {
				if (matching(i, j) == 1) {
					logp += likelihoods(i, j);
					matched_durations(i) = simulated_durations(used_durations[j]);
				}
			}
		}
	}

	return true;
}

void DataLikelihoodDuration::NotifyDuration(size_t cell_ix, Real duration)
{
	simulated_durations(cell_ix) = duration;
}
