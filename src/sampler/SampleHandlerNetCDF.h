#pragma once

#include "NetCDFDataFile.h"
#include "SamplerPT.h"

namespace bcm3 {

class SampleHandlerNetCDF : public SampleHandler
{
public:
	SampleHandlerNetCDF();
	~SampleHandlerNetCDF();
	
	void SetFile(const std::string& filename);
	bool Initialize(size_t sample_count, const VariableSet* varset, VectorReal output_temperatures);
	void Close();

	virtual void ReceiveSample(const VectorReal& values, Real lprior, Real llh, Real temperature);

private:
	std::string filename;
	NetCDFDataFile netcdf_file;
	VectorReal output_temperatures;
	size_t num_variables;
	std::vector<size_t> sample_ix;
	bool keep_output_open;
};

}
