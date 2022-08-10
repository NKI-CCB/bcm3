#pragma once

#include "SamplerPT.h"

namespace bcm3 {

class SampleHandlerTSV : public SampleHandler
{
public:
	SampleHandlerTSV();
	~SampleHandlerTSV();
	
	void SetFile(const std::string& filename);
	bool Initialize(size_t sample_count, const std::vector<std::string>& variables, VectorReal output_temperatures);

	virtual void ReceiveSample(const VectorReal& values, Real lprior, Real llh, Real temperature);

private:
	std::string filename;
	size_t num_variables;
};

}
