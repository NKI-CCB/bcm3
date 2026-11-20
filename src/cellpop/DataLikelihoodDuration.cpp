#include "DataLikelihoodDuration.h"
#include "Experiment.h"
#include "HungarianMatching.h"
#include "ProbabilityDistributions.h"

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
	experiment->AddRequestedDuration(this, period);
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
}

bool DataLikelihoodDuration::Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	Real stdev = GetCurrentSTDev(transformed_values, non_sampled_parameters);
	logp = 0.0;

	std::vector<Eigen::Index> observed_indices, simulated_indices;
	for (Eigen::Index i = 0; i < observed_durations.size(); i++) {
		if (std::isfinite(observed_durations(i))) {
			observed_indices.push_back(i);
		}
	}
	for (Eigen::Index j = 0; j < simulated_durations.size(); j++) {
		if (std::isfinite(simulated_durations(j))) {
			simulated_indices.push_back(j);
		}
	}
	if (observed_indices.empty()) {
		return true;
	}
	if (observed_indices.size() > simulated_indices.size()) {
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	}

	MatrixReal likelihoods(observed_indices.size(), simulated_indices.size());
	for (size_t i = 0; i < observed_indices.size(); i++) {
		for (size_t j = 0; j < simulated_indices.size(); j++) {
			likelihoods(i, j) = bcm3::LogPdfNormal(observed_durations(observed_indices[i]), simulated_durations(simulated_indices[j]), stdev);
		}
	}

	auto matching = bcm3::RectangularHungarian(-likelihoods);
	if (!matching) {
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	}

	for (size_t i = 0; i < observed_indices.size(); i++) {
		int j = (*matching)[i];
		if (j == -1) {
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}
		logp += likelihoods(i, j);
		matched_durations(observed_indices[i]) = simulated_durations(simulated_indices[j]);
	}
	logp *= weight;

	return true;
}

void DataLikelihoodDuration::NotifyDuration(size_t cell_ix, Real duration)
{
	simulated_durations(cell_ix) = duration;
}
