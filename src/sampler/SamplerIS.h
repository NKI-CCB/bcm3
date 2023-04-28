#pragma once

#include "Sampler.h"

namespace bcm3 {

	class SamplerIS : public Sampler
	{
	public:
		SamplerIS(size_t threads, size_t max_memory_use);
		~SamplerIS();

		virtual bool LoadSettings(const boost::program_options::variables_map& vm);
		virtual bool Initialize();

		static void AddOptionsDescription(boost::program_options::options_description& pod);

	private:
		virtual bool RunImpl();
		virtual void LogStatistics();

		Real heighest_weight;
	};

}
