#include "Utils.h"
#include "LikelihoodMitosisTimeEstimation.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"
//#include "../../dependencies/HungarianAlgorithm-master/matching.h"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/sobol.hpp>

LikelihoodMitosisTimeEstimation::LikelihoodMitosisTimeEstimation(size_t numthreads, size_t evaluation_threads)
{
}

LikelihoodMitosisTimeEstimation::~LikelihoodMitosisTimeEstimation()
{
}

bool LikelihoodMitosisTimeEstimation::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;

	bcm3::NetCDFDataFile data;
	std::string data_filename = likelihood_node.get<std::string>("<xmlattr>.data_file", "trajectories.nc");
	if (!data.Open(data_filename, false)) {
		return false;
	}

	if (!data.GetMatrix("simulation", "trajectories", observed_trajectories)) {
		return false;
	}
	const size_t ncell = observed_trajectories.cols();
	if (ncell == 0) {
		LOGERROR("No cells");
		return false;
	}

	if (!data.GetVector("simulation", "time", timepoints)) {
		return false;
	}
	if (timepoints.size() == 0) {
		LOGERROR("No timepoints");
		return false;
	}

	boost::random::uniform_01<Real> distribution;
	boost::random::sobol sobol_sequence(2);

	sobol_sequence_values.resize(ncell, 2);
	for (size_t i = 0; i < ncell; i++) {
		for (int j = 0; j < 2; j++) {
			Real val = distribution(sobol_sequence);
			sobol_sequence_values(i, j) = pow(2.0, bcm3::QuantileNormal(val, 0, 0.5));
		}
	}

	simulated_trajectories.resize(timepoints.size(), ncell);
	likelihoods.resize(ncell, ncell);

	return true;
}

bool LikelihoodMitosisTimeEstimation::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	threadix = 0;

	const Real mitosis_times_stdev = bcm3::fastpow10(values[varset->GetVariableIndex("mitosis_times_stdev")]);
	const Real entry_time_stdev = bcm3::fastpow10(values[varset->GetVariableIndex("entry_time_stdev")]);
	const Real trajectory_noise_stdev = bcm3::fastpow10(values[varset->GetVariableIndex("trajectory_noise_stdev")]);
	const size_t ncell = observed_trajectories.cols();

	VectorReal simulated_times = sobol_sequence_values.col(0) * mitosis_times_stdev;
	VectorReal start_times = sobol_sequence_values.col(1) * entry_time_stdev;

	// Generate trajectories
	simulated_trajectories.setZero();
	for (size_t i = 0; i < ncell; i++) {
		for (size_t j = 0; j < timepoints.size(); j++) {
			if (timepoints(j) >= start_times(i) && timepoints(j) < start_times(i) + simulated_times(i)) {
				simulated_trajectories(j, i) = 1.0;
			}
		}
	}

	Real inv_two_sigma = 1.0 / (2.0 * trajectory_noise_stdev * trajectory_noise_stdev);
	Real C = -log(trajectory_noise_stdev) - 0.91893853320467274178032973640562;
	for (size_t i = 0; i < ncell; i++) {
		for (size_t j = 0; j < ncell; j++) {
			Real match_logp = 0.0;
			for (size_t k = 0; k < timepoints.size(); k++) {
				Real d = observed_trajectories(k, j) - simulated_trajectories(k, i);
				match_logp -= d * d * inv_two_sigma;
			}
			match_logp += timepoints.size() * C;
			likelihoods(i, j) = match_logp;
		}
	}

#if TODO
	MatrixReal matching;
	if (!HungarianMatching(likelihoods, matching, HUNGARIAN_MATCH_MAX)) {
		LOGERROR("Could not find matching");
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	}

	logp = 0.0;
	for (size_t i = 0; i < ncell; i++) {
		for (size_t j = 0; j < ncell; j++) {
			if (matching(i, j) == 1) {
				logp += likelihoods(i, j);
			}
		}
	}
#endif

	return true;
}
