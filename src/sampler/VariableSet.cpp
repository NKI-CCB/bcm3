#include "Utils.h"
#include "VariableSet.h"

#include <boost/property_tree/xml_parser.hpp>

namespace bcm3 {

VariableSet::VariableSet()
{
}

VariableSet::~VariableSet()
{
}

bool VariableSet::LoadFromXML(const std::string& filename)
{
	Name = filename;
	Name.at(Name.find_last_of('.')) = '_';

	boost::property_tree::ptree pt;
	try {
		boost::property_tree::read_xml(filename, pt);
	} catch (boost::property_tree::xml_parser_error &e) {
		LOGERROR("Error loading variable file: %s", e.what());
		return false;
	}
	
	try {
		boost::property_tree::ptree prior_node;
		if (pt.count("prior") > 0) {
			prior_node = pt.get_child("prior");
		} else if (pt.count("variableset") > 0) {
			prior_node = pt.get_child("variableset");
		} else {
			LOGERROR("Incorrect prior XML format");
			return false;
		}
		
		BOOST_FOREACH(const boost::property_tree::ptree::value_type& var, prior_node.get_child("")) {
			if (var.first == "variable") {
				std::string name = var.second.get<std::string>("<xmlattr>.name");
				variables.push_back(name);

				bool logspace = var.second.get<bool>("<xmlattr>.logspace", false);
				bool logistic = var.second.get<bool>("<xmlattr>.logistic", false);

				if (logspace) {
					transforms.push_back(Transform_Log10);
				} else if (logistic) {
					transforms.push_back(Transform_Logit);
				} else {
					transforms.push_back(Transform_None);
				}
			}
		}

		return true;
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing variable file: %s", e.what());
		return false;
	}
}

void VariableSet::AddVariable(const std::string& name, bool logspace, bool logistic)
{
	variables.push_back(name);

	if (logspace) {
		transforms.push_back(Transform_Log10);
	} else if (logistic) {
		transforms.push_back(Transform_Logit);
	} else {
		transforms.push_back(Transform_None);
	}
}

size_t VariableSet::GetVariableIndex(const std::string& name, bool log_error) const
{
	for (size_t vi = 0; vi < variables.size(); vi++) {
		if (variables[vi] == name) {
			return vi;
		}
	}
	if (log_error) {
		LOGERROR("Could not find variable \"%s\"", name.c_str());
	}
	return std::numeric_limits<size_t>::max();
}

Real VariableSet::TransformVariable(size_t i, Real x) const
{
	ASSERT(i < variables.size());
	ASSERT(i < transforms.size());

	switch (transforms[i]) {
	case Transform_None:
		return x;

	case Transform_Log:
		return exp(x);

	case Transform_Log10:
		return bcm3::fastpow10(x);

	case Transform_Logit:
		if (x > 0) {
			Real z = exp(-x);
			return 1.0 / (1.0 + z);
		} else {
			Real z = exp(x);
			return z / (1.0 + z);
		}

	default:
		return x;
	}
}

}
