#pragma once

class Experiment;

class TreatmentTrajectory
{
public:
	bool Load(Experiment* experiment, const bcm3::NetCDFDataFile& data_file, const std::string& treatment_variable);
	Real GetConcentration(Real time, size_t sample_ix);
	size_t GetNumSamples() { return concentrations.rows(); }

private:
	VectorReal timepoints;
	MatrixReal concentrations;
};
