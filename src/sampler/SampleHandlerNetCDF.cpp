#include "Utils.h"
#include "NetCDFDataFile.h"
#include "SampleHandlerNetCDF.h"
#include "VariableSet.h"

namespace bcm3 {

	SampleHandlerNetCDF::SampleHandlerNetCDF()
		: num_variables(0)
		, sample_ix(0)
		, keep_output_open(true)
	{
	}

	SampleHandlerNetCDF::~SampleHandlerNetCDF()
	{
	}

	void SampleHandlerNetCDF::SetFile(const std::string& filename)
	{
		this->filename = filename;
	}

	bool SampleHandlerNetCDF::Initialize(size_t sample_count, const VariableSet* varset, VectorReal output_temperatures)
	{
		if (!netcdf_file.Create(filename)) {
			LOGERROR("Cannot open output file");
			return false;
		}

		num_variables = varset->GetNumVariables();

		this->output_temperatures = output_temperatures;
		std::vector<unsigned int> sampleix(sample_count);
		for (size_t i = 0; i < sample_count; i++) {
			sampleix[i] = (unsigned int)i+1;
		}
		bool result = true;

		std::vector<std::string> varnames(num_variables);
		for (size_t vi = 0; vi < varset->GetNumVariables(); vi++) {
			varnames[vi] = varset->GetVariableName(vi);
		}

		// Allocate variables
		result &= netcdf_file.CreateGroup("samples");
		result &= netcdf_file.CreateDimension("samples", "sample_ix", sampleix);
		result &= netcdf_file.CreateDimension("samples", "variable", varnames);
		result &= netcdf_file.CreateDimension("samples", "temperature", output_temperatures);
		result &= netcdf_file.CreateVariable<unsigned int>("samples", "variable_transform", "variable");
		result &= netcdf_file.CreateVariable("samples", "variable_values", "sample_ix", "temperature", "variable");
		result &= netcdf_file.CreateVariable("samples", "log_prior", "sample_ix", "temperature");
		result &= netcdf_file.CreateVariable("samples", "log_likelihood", "sample_ix", "temperature");
		result &= netcdf_file.CreateVariable("samples", "weights", "sample_ix", "temperature");

		for (size_t vi = 0; vi < varset->GetNumVariables(); vi++) {
			result &= netcdf_file.PutValue("samples", "variable_transform", vi, (unsigned int)varset->GetVariableTransform(vi));
		}

		if (!keep_output_open) {
			netcdf_file.Close();
		}

		sample_ix.resize(output_temperatures.size(), 0);
		return result;
	}

	void SampleHandlerNetCDF::Close()
	{
		netcdf_file.Close();
	}

	void SampleHandlerNetCDF::ReceiveSample(const VectorReal& values, Real lprior, Real llh, Real temperature, Real weight)
	{
		if (!keep_output_open) {
			if (!netcdf_file.Open(filename, true)) {
				LOGERROR("Cannot open output file");
				return;
			}
		}

		size_t temperature_ix = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < output_temperatures.size(); i++) {
			if (temperature == output_temperatures[i]) {
				temperature_ix = i;
				break;
			}
		}
		if (temperature_ix != std::numeric_limits<size_t>::max()) {
			size_t si = sample_ix[temperature_ix];

			netcdf_file.PutValue("samples", "sample_ix", si, (unsigned int)sample_ix[temperature_ix]);
			netcdf_file.PutValue("samples", "log_prior", si, temperature_ix, lprior);
			netcdf_file.PutValue("samples", "log_likelihood", si, temperature_ix, llh);
			netcdf_file.PutValue("samples", "weights", si, temperature_ix, weight);
			netcdf_file.PutValuesDim3("samples", "variable_values", si, temperature_ix, 0, values);

			if (!keep_output_open) {
				netcdf_file.Close();
			}

			sample_ix[temperature_ix]++;

			if (temperature_ix == sample_ix.size() - 1 && si % 10 == 0) {
				netcdf_file.Sync();
			}
		} else {
			// This temperature is not stored in the output value
		}
	}

}
