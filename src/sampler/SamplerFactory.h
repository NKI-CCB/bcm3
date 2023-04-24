#pragma once

#include "Sampler.h"
#include <boost/program_options.hpp>

namespace bcm3 {

	class Sampler;

	class SamplerFactory
	{
	public:
		static std::shared_ptr<bcm3::Sampler> Create(const boost::program_options::variables_map& vm, size_t numthreads, size_t max_memory_use);
		static void AddOptionsDescription(boost::program_options::options_description& pod);
	};

}
