#include "Utils.h"
#include "Prior.h"
#include "PriorIndependence.h"
#include "VariableSet.h"

#include <boost/property_tree/xml_parser.hpp>

namespace bcm3 {

	Prior::Prior(std::shared_ptr<VariableSet> varset)
		: varset(varset)
	{
		LowerBounds = VectorReal::Constant(varset->GetNumVariables(), -std::numeric_limits<Real>::infinity());
		UpperBounds = VectorReal::Constant(varset->GetNumVariables(), std::numeric_limits<Real>::infinity());
	}

	Prior::~Prior()
	{
	}

	std::unique_ptr<Prior> Prior::Create(const std::string& filename, std::shared_ptr<VariableSet> varset, size_t numthreads)
	{
		std::unique_ptr<Prior> prior;

		boost::property_tree::ptree pt;
		try {
			boost::property_tree::read_xml(filename, pt);
		} catch (boost::property_tree::xml_parser_error &e) {
			LOGERROR("Error loading variable file: %s", e.what());
			return prior;
		}

		try {
			boost::property_tree::ptree prior_node;
			if (pt.count("prior") > 0) {
				prior_node = pt.get_child("prior");
			} else if (pt.count("variableset") > 0) {
				prior_node = pt.get_child("variableset");
			} else {
				LOGERROR("Incorrect prior XML format");
				return prior;
			}

			std::string type = prior_node.get<std::string>("<xmlattr>.type", "");

			if (type.empty() || type == "independence") {
				prior = std::make_unique<PriorIndependence>(varset);
			} else {
				LOGERROR("Unknown prior type \"%s\"", type.c_str());
				return prior;
			}

			if (prior->LoadFromXML(prior_node)) {
				return prior;
			} else {
				prior.reset();
			}
		} catch (boost::property_tree::ptree_error &e) {
			LOGERROR("Error parsing variable file: %s", e.what());
			prior.reset();
		}

		return prior;
	}

}
