#include "Utils.h"
#include "SampleHandlerTSV.h"

namespace bcm3 {

SampleHandlerTSV::SampleHandlerTSV()
	: num_variables(0)
{
}

SampleHandlerTSV::~SampleHandlerTSV()
{
}

void SampleHandlerTSV::SetFile(const std::string& filename)
{
	this->filename = filename;
}

bool SampleHandlerTSV::Initialize(size_t sample_count, const std::vector<std::string>& variables, VectorReal output_temperatures)
{
	FILE* file = fopen(filename.c_str(), "w");
	if (!file) {
		LOGERROR("Failed to open output file \"%s\"", filename.c_str());
		return false;
	} else {
		for (auto vi : variables) {
			fprintf(file, "%s\t", vi.c_str());
		}
		fprintf(file, "log prior\tlog likelihood\tweight\n");
		fclose(file);
	}
	return true;
}


void SampleHandlerTSV::ReceiveSample(const VectorReal& values, Real lprior, Real llh, Real temperature, Real weight)
{
	if (temperature == 1.0) {
		FILE* file = fopen(filename.c_str(), "a");
		if (!file) {
			LOGERROR("Failed to open output file \"%s\"", filename.c_str());
		} else {
			for (size_t i = 0; i < values.size(); i++) {
				fprintf(file, "%.6g\t", values(i));
			}
			fprintf(file, "%.6g\t", lprior);
			fprintf(file, "%.6g\n", llh);
			fprintf(file, "%.6g\n", weight);
			fclose(file);
		}
	}
}

}
