#include "Utils.h"
#include "CSVParser.h"
#include "LikelihoodCellCycleMarker.h"
#include "ProbabilityDistributions.h"

//const Real global_data[] = {
//	6.24,6.23,6.27,6.11,6.28,6.14,6.17,6.13,6,5.97,6.1,6.07,6.02,6.1,6.03,6.05,6.08,5.95,6.07,6.13,6.12,6.23,6.16,6.19,6.13,5.9,6.19,6.2,6.24,6.37,6.46,6.47,6.38,6.48,6.71,6.97,7.46,8.02,8.59,9.51,10.52,11.59,12.52,13.6,14.76,16.53,17.53,19.48,20.98,21.95,23.47,24.62,27.09,28.96,30.44,32.55,33.29,34.42,36.56,38.26,40.19,41.63,43.82,45.79,48.03,48.67,51.08,52.27,53.45,55.13,56.6,58.08,59.17,59.47,61.07,60.98,62.64,63.58,64.13,64.68,66.51,67.45,68.83,69.18,70.59,71.26,71.5,72.19,74.19,74.75,74.47,73.82,76.91,77.02,78.62,80,80.8,80.77,82.22,82.23,83.32,82.98,84,83.6,85.36,86.65,86.03,85.59,87.17,86.6,88.32,88.2,88.23,88.65,89.11,88.95,88.39,88.45,88.67,87.95,87.36,87.96,85.93,85.98,85.77,86.89,86.68,88.56,89.24,88.84,89.43,87.82,87.91,87.56,89.18,87.25,87.85,85.91,85.85,87.35,85.61,85.99,88.65,89.48,88.65,86.74,85.04,85.21,85.87,86.76,83.59,82.94,84.83,85.28,86.63,89.45,90.39,88.56,91.08,91.2,85.76,89.13,91.27,90.57,85.89,89.23,86.31,80.89,80.36,83.28,79.32,79.68,78.27,81,74.85,53.38,31.48,32.73,32.32,30.47,29.3,27.77,26.29,25.39,23.32,23,23.54,23.46,23.21,22.95,22.29,22.34,21.45,20.2,19.39,19.29,19.97,19.23,18.93,18.99,18.33,18.02,18.05,17.84,45.09,std::numeric_limits<Real>::quiet_NaN(),18.42,std::numeric_limits<Real>::quiet_NaN(),17.24,16.81,17.16,17.1,16.25,15.74,16.06,16.11,15.97,15.55,15.13,14.42,13.97,13.68,13.24,13.32
//};

LikelihoodCellCycleMarker::LikelihoodCellCycleMarker(size_t sampling_threads, size_t evaluation_threads)
{
}

LikelihoodCellCycleMarker::~LikelihoodCellCycleMarker()
{
}

bool LikelihoodCellCycleMarker::Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm)
{
	this->varset = varset;
	if (varset->GetNumVariables() != 10) {
		LOGERROR("Variable set should contain exactly 10 variables but contains %zd variables.", varset->GetNumVariables());
		return false;
	}

	int use_track_ix = vm["ccm.track_ix"].as<int>();

	try {
		std::string data_file = likelihood_node.get<std::string>("<xmlattr>.data_file");

		bcm3::CSVParser parser;
		if (!parser.Parse(data_file, "\t")) {
			LOGERROR("Failed to parse data file \"%s\"", data_file.c_str());
			return false;
		}

		if (use_track_ix < 0 || use_track_ix >= parser.GetNumRows()) {
			LOGERROR("Track index %d is out of bounds for data file with %zd columns.", use_track_ix, parser.GetNumColumns());
			return false;
		} else {
			LOG("Using track index %d.", use_track_ix);
		}

		data = VectorReal(parser.GetNumColumns());
		for (size_t i = 0; i < parser.GetNumColumns(); i++) {
			data(i) = parser.GetEntry(use_track_ix, i);
		}
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return true;
}

bool LikelihoodCellCycleMarker::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp)
{
	Real S_entry_time = values(0);
	Real S_duration = values(1);
	Real plateau_duration = values(2);
	Real plateau_time = S_entry_time + S_duration;
	Real mitosis_time = S_entry_time + S_duration + plateau_duration;
	Real base_signal = values(3);
	Real S_signal_increase = values(4);
	Real plateau_signal_increase = values(5);
	Real mitosis_signal_fraction = values(6);
	Real mitosis_signal_decrease = values(7);
	Real additive_noise = values(8);
	Real proportional_noise = values(9);

	logp = 0.0;
	for (int i = 0; i < data.size(); i++) {
		Real x = base_signal;
		if (i > mitosis_time) {
			x += (S_duration * S_signal_increase + plateau_duration * plateau_signal_increase) * mitosis_signal_fraction;
			x -= mitosis_signal_decrease * (i - mitosis_time);
		} else if (i > plateau_time) {
			x += S_duration * S_signal_increase + (i - plateau_time) * plateau_signal_increase;
		} else if (i > S_entry_time) {
			x += S_signal_increase * (i - S_entry_time);
		}

		Real y = data[i];

		Real sigma = additive_noise + proportional_noise * std::max(x, 0.0);
		logp += bcm3::LogPdfTnu4(y, x, sigma, true);
	}

	return true;
}

void LikelihoodCellCycleMarker::AddOptionsDescription(boost::program_options::options_description& pod)
{
	pod.add_options()
		("ccm.track_ix", boost::program_options::value<int>()->default_value(0), "Use the specified trace in the data file.")
	;
}
