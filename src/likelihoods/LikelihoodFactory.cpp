#include "Utils.h"
#include "LikelihoodFactory.h"

#include "LikelihoodDLL.h"
#include "LikelihoodDummy.h"
#include "LikelihoodIncucytePopulation.h"
#include "LikelihoodMitosisTimeEstimation.h"
#include "LikelihoodODE.h"
#include "LikelihoodPharmacokineticTrajectory.h"
#include "LikelihoodPopPKTrajectory.h"
#include "PharmacoLikelihoodPopulation.h"
#include "PharmacoLikelihoodSingle.h"
#include "TestLikelihoodBanana.h"
#include "TestLikelihoodCircular.h"
#include "TestLikelihoodMultimodalGaussians.h"
#include "TestLikelihoodTruncatedT.h"

#if BCM_INCLUDE_fISA
#include "fISALikelihood.h"
#endif

#if BCM_INCLUDE_CELLPOP
#include "CellPopulationLikelihood.h"
#endif

#include <boost/property_tree/xml_parser.hpp>

namespace bcm3 {

std::shared_ptr<bcm3::Likelihood> LikelihoodFactory::CreateLikelihood(std::string likelihood_xml_fn, std::shared_ptr<const bcm3::VariableSet> varset, const boost::program_options::variables_map& vm, size_t sampling_threads, size_t evaluation_threads, bool running_inference)
{
	std::shared_ptr<bcm3::Likelihood> ll;

	// Load likelihood file
	boost::property_tree::ptree pt;
	try {
		boost::property_tree::read_xml(likelihood_xml_fn, pt);
	} catch (boost::property_tree::xml_parser_error & e) {
		LOGERROR("Error loading likelihood file: %s", e.what());
		return ll;
	}

	try {
		boost::property_tree::ptree likelihood_node = pt.get_child("bcm_likelihood");
		
		std::string type = likelihood_node.get<std::string>("<xmlattr>.type");
		if (type == "dll") {
			ll = std::make_shared<LikelihoodDLL>(sampling_threads, evaluation_threads);
		} else if (type == "dummy") {
			ll = std::make_shared<LikelihoodDummy>(sampling_threads, evaluation_threads);
		} else if (type == "incucyte_population") {
			ll = std::make_shared<LikelihoodIncucytePopulation>(sampling_threads, evaluation_threads);
		} else if (type == "mitosis_time_estimation") {
			ll = std::make_shared< LikelihoodMitosisTimeEstimation>(sampling_threads, evaluation_threads);
		} else if (type == "ODE") {
			ll = std::make_shared<LikelihoodODE>(sampling_threads, evaluation_threads);
		} else if (type == "pharmacokinetic_trajectory") {
			ll = std::make_shared<LikelihoodPharmacokineticTrajectory>(sampling_threads, evaluation_threads);
		} else if (type == "pop_pk_trajectory") {
			ll = std::make_shared<LikelihoodPopPKTrajectory>(sampling_threads, evaluation_threads);
		} else if (type == "pharmaco_single") {
			ll = std::make_shared<PharmacoLikelihoodSingle>(sampling_threads, evaluation_threads);
		} else if (type == "pharmaco_population") {
			ll = std::make_shared<PharmacoLikelihoodPopulation>(sampling_threads, evaluation_threads);
		} else if (type == "banana") {
			ll = std::make_shared<TestLikelihoodBanana>(sampling_threads, evaluation_threads);
		} else if (type == "circular") {
			ll = std::make_shared<TestLikelihoodCircular>(sampling_threads, evaluation_threads);
		} else if (type == "multimodal_gaussians") {
			ll = std::make_shared<TestLikelihoodMultimodalGaussians>(sampling_threads, evaluation_threads);
		} else if (type == "truncated_t") {
			ll = std::make_shared<TestLikelihoodTruncatedT>(sampling_threads, evaluation_threads);
#if INCLUDE_fISA
		} else if (type == "fISA") {
			ll = std::make_shared<fISALikelihood>(sampling_threads, evaluation_threads);
#endif
#if BCM_INCLUDE_CELLPOP
		} else if (type == "cell_population") {
			ll = std::make_shared<CellPopulationLikelihood>(sampling_threads, evaluation_threads, !running_inference);
#endif
		} else {
			LOGERROR("Unknown likelihood type \"%s\"", type.c_str());
		}
		if (ll) {
			if (!ll->Initialize(varset, likelihood_node, vm)) {
				LOGERROR("Failed to initialize likelihood");
				ll.reset();
			}
		} else {
			LOGERROR("Failed to create likelihood");
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return ll;
	}

	return ll;
}

void LikelihoodFactory::AddOptionsDescription(boost::program_options::options_description& pod)
{
#if BCM_INCLUDE_CELLPOP
	CellPopulationLikelihood::AddOptionsDescription(pod);
#endif
	LikelihoodPharmacokineticTrajectory::AddOptionsDescription(pod);
	PharmacoLikelihoodSingle::AddOptionsDescription(pod);
}

}
