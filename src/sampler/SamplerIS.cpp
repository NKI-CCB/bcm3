#include "Utils.h"
#include "SamplerIS.h"
#include "VariableSet.h"

namespace bcm3 {

	SamplerIS::SamplerIS(size_t threads, size_t max_memory_use)
		: Sampler(threads, max_memory_use)
	{
	}

	SamplerIS::~SamplerIS()
	{
	}

	bool SamplerIS::LoadSettings(const boost::program_options::variables_map& vm)
	{
		if (!Sampler::LoadSettings(vm)) {
			return false;
		}

		try {
		} catch (boost::program_options::error& e) {
			LOGERROR("Error parsing sampler: %s", e.what());
			return false;
		}

		temperatures.setConstant(1, 1.0);

		return true;
	}

	bool SamplerIS::Initialize()
	{
		if (!Sampler::Initialize()) {
			return false;
		}

		return true;
	}

	void SamplerIS::AddOptionsDescription(boost::program_options::options_description& pod)
	{
	}

	bool SamplerIS::RunImpl()
	{
		bool result = true;

		LOG("Starting sampling loop...");
		UpdateProgress(0.0, true);
		size_t total_samples = num_samples * use_every_nth;
		for (size_t si = 0; si < total_samples; si++) {
			UpdateProgress(si / (Real)total_samples, false);

			VectorReal sample;
			if (!prior->Sample(sample, &rng)) {
				return false;
			}

			Real lprior, llh;
			if (!EvaluatePriorLikelihood(0, sample, lprior, llh)) {
				return false;
			}

			Real lposterior = lprior + llh;
			Real lweight = lposterior - lprior;
			if (lweight > -10) {
				int a;
				a = 6;
			}
			for (auto handler : sample_handlers) {
				handler->ReceiveSample(sample, lprior, llh, 1.0, exp(lweight));
			}
		}

		return true;
	}

	void SamplerIS::LogStatistics()
	{
	}

}
