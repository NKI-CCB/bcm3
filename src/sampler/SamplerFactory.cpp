#include "Utils.h"
#include "SamplerFactory.h"
#include "SamplerIS.h"
#include "SamplerPT.h"

namespace po = boost::program_options;

namespace bcm3 {

	std::shared_ptr<bcm3::Sampler> SamplerFactory::Create(const boost::program_options::variables_map& vm, size_t threads, size_t max_memory_use)
	{
		std::shared_ptr<bcm3::Sampler> sampler;

		std::string type;
		try {
			type = vm["sampler.type"].as<std::string>();
		} catch (boost::program_options::error& e) {
			LOGERROR("Error parsing sampler config: %s", e.what());
			return NULL;
		}

		if (type == "ptmh" || type == "parallel_tempered_Metropolis_Hastings") {
			sampler = std::make_shared<SamplerPT>(threads, max_memory_use);
		} else if (type == "is" || type == "importance_sampling") {
			sampler = std::make_shared<SamplerIS>(threads, max_memory_use);
		}

		if (!sampler) {
			return sampler;
		}

		if (sampler->LoadSettings(vm)) {
			return sampler;
		} else {
			sampler.reset();
			return sampler;
		}
	}

	void SamplerFactory::AddOptionsDescription(boost::program_options::options_description& pod)
	{
		pod.add_options()
			("sampler.type", boost::program_options::value<std::string>()->default_value("ptmh"), "Sampling algorithm to use [ptmh|is].")
			;

		Sampler::AddOptionsDescription(pod);
		SamplerIS::AddOptionsDescription(pod);
		SamplerPT::AddOptionsDescription(pod);
	}
}
