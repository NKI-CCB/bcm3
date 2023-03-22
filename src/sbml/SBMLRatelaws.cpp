#include "Utils.h"
#include "SBMLRatelaws.h"

#define USE_INTEGER_HILL_OPTIMIZATIONS 1

inline Real hill_function(Real x, Real k, Real n)
{
	Real xn = pow(x, n);
	Real kn = pow(k, n);
	return xn / (kn + xn);
}

inline Real michaelis_menten_function(Real kcat, Real KM, Real e, Real s)
{
	if (s < 0) return -kcat * e * s / (-KM + s);
	if (s < DBL_EPSILON) return 0.0;
	return kcat * e * s / (KM + s);
}

inline Real michaelis_menten_derivative_enzyme(Real kcat, Real KM, Real e, Real s)
{
	if (s < 0) return -kcat * s / (-KM + s);
	if (s < DBL_EPSILON) return 0.0;
	return kcat * s / (KM + s);
}

inline Real michaelis_menten_derivative_substrate(Real kcat, Real KM, Real e, Real s)
{
	if (s < 0) return e * kcat * KM / (bcm3::square(-KM + s));
	if (s < DBL_EPSILON) return 0.0;
	return e * kcat * KM / (bcm3::square(KM + s));
}

inline Real safepow(Real x, Real n)
{
	if (x < 0) {
		return 0.0;
	} else {
		return pow(x, n);
	}
}

bool SBMLRatelawElement::Generate(const ASTNode* node,
								  const Model* model,
								  const std::map<std::string, size_t>& species_index_map,
								  const std::map<std::string, size_t>& parameter_index_map,
								  const std::map<std::string, size_t>& constant_species_index_map,
								  const std::map<std::string, size_t>& non_sampled_parameter_index_map,
								  const std::map<std::string, Real>& forced_parameter_values,
								  std::unique_ptr<SBMLRatelawElement>* element)
{
	switch (node->getType()) {
	case AST_PLUS:
		*element = std::make_unique<SBMLRatelawElementPlus>();
		(*element)->children.resize(node->getNumChildren());
		for (unsigned int i = 0; i < node->getNumChildren(); i++) {
			if (!Generate(node->getChild(i), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[i])) {
				return false;
			}
		}
		break;

	case AST_MINUS:
		if (node->getNumChildren() == 1) {
			*element = std::make_unique<SBMLRatelawElementNegate>();
			(*element)->children.resize(1);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
		}
		else if (node->getNumChildren() == 2) {
			*element = std::make_unique<SBMLRatelawElementMinus>();
			(*element)->children.resize(2);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
			if (!Generate(node->getChild(1), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[1])) {
				return false;
			}
		}
		else {
			// I'm not sure how this is defined..
			LOGERROR("Not implemented");
			return false;
		}
		break;

	case AST_TIMES:
		*element = std::make_unique<SBMLRatelawElementTimes>();
		(*element)->children.resize(node->getNumChildren());
		for (unsigned int i = 0; i < node->getNumChildren(); i++) {
			if (!Generate(node->getChild(i), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[i])) {
				return false;
			}
		}
		break;

	case AST_DIVIDE:
		if (node->getNumChildren() == 2) {
			*element = std::make_unique<SBMLRatelawElementDivide>();
			(*element)->children.resize(2);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
			if (!Generate(node->getChild(1), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[1])) {
				return false;
			}
		}
		else {
			// I'm not sure how this is defined..
			LOGERROR("AST division with more than 2 children?");
			return false;
		}
		break;

	case AST_NAME:
		{
			const char* name = node->getName();
			std::map<std::string, size_t>::const_iterator iter;

			// Look for a forced parameter value
			std::map<std::string, Real>::const_iterator fiter = forced_parameter_values.find(name);
			if (fiter != forced_parameter_values.end()) {
				std::unique_ptr< SBMLRatelawElementConstant> el = std::make_unique<SBMLRatelawElementConstant>();
				el.get()->constant = fiter->second;
				*element = std::move(el);
				break;
			}

			// Look for an inference parameter
			iter = parameter_index_map.find(name);
			if (iter != parameter_index_map.end()) {
				std::unique_ptr< SBMLRatelawElementLookupParameter> el = std::make_unique<SBMLRatelawElementLookupParameter>();
				el.get()->ix = iter->second;
				*element = std::move(el);
				break;
			}

			// Look for a species
			iter = species_index_map.find(name);
			if (iter != species_index_map.end()) {
				std::unique_ptr< SBMLRatelawElementLookupSpecies> el = std::make_unique<SBMLRatelawElementLookupSpecies>();
				el.get()->ix = iter->second;
				*element = std::move(el);
				break;
			}

			// Look for a constant species
			iter = constant_species_index_map.find(name);
			if (iter != constant_species_index_map.end()) {
				std::unique_ptr< SBMLRatelawElementLookupConstantSpecies> el = std::make_unique<SBMLRatelawElementLookupConstantSpecies>();
				el.get()->ix = iter->second;
				*element = std::move(el);
				break;
			}

			// Look for a non-sampled parameter
			iter = non_sampled_parameter_index_map.find(name);
			if (iter != non_sampled_parameter_index_map.end()) {
				std::unique_ptr< SBMLRatelawElementLookupNonSampledParameter> el = std::make_unique<SBMLRatelawElementLookupNonSampledParameter>();
				el.get()->ix = iter->second;
				*element = std::move(el);
				break;
			}

			// If there wasn't an inference parameter, look if the value has been specified in the SBML document
			bool found_parameter = false;
			for (unsigned int i = 0; i < model->getNumParameters(); i++) {
				if (model->getParameter(i)->getId() == name) {
					std::unique_ptr< SBMLRatelawElementConstant> el = std::make_unique<SBMLRatelawElementConstant>();
					el.get()->constant = model->getParameter(i)->getValue();
					*element = std::move(el);
					found_parameter = true;
					break;
				}
			}
			if (found_parameter) {
				break;
			}

			LOGERROR("AST_NAME name \"%s\" does not map to either a species id or a parameter", name);
			return false;
		}
		break;

	case AST_INTEGER:
		{
			std::unique_ptr< SBMLRatelawElementConstant> el = std::make_unique<SBMLRatelawElementConstant>();
			el.get()->constant = node->getInteger();
			*element = std::move(el);
		}
		break;

	case AST_REAL:
	case AST_REAL_E:
		{
			std::unique_ptr< SBMLRatelawElementConstant> el = std::make_unique<SBMLRatelawElementConstant>();
			el.get()->constant = node->getReal();
			*element = std::move(el);
		}
		break;

	case AST_FUNCTION_EXP:
		if (node->getNumChildren() == 1) {
			*element = std::make_unique<SBMLRatelawElementFunctionExp>();
			(*element)->children.resize(1);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
		} else {
			LOGERROR("AST exp with wrong number of children?");
			return false;
		}
		break;

	case AST_FUNCTION_LN:
		if (node->getNumChildren() == 1) {
			*element = std::make_unique<SBMLRatelawElementFunctionLog>();
			(*element)->children.resize(1);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
		} else {
			LOGERROR("AST log with wrong number of children?");
			return false;
		}
		break;

	case AST_FUNCTION_POWER:
		if (node->getNumChildren() == 2) {
			*element = std::make_unique<SBMLRatelawElementFunctionPow>();
			(*element)->children.resize(2);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
			if (!Generate(node->getChild(1), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[1])) {
				return false;
			}
		} else {
			LOGERROR("AST pow with wrong number of children?");
			return false;
		}
		break;

	case AST_FUNCTION:
		if (strcmp(node->getName(), "hill") == 0) {
			if (node->getNumChildren() != 3) {
				LOGERROR("Hill function should have three parameters");
				return false;
			}
			*element = std::make_unique<SBMLRatelawElementFunctionHill>();
			(*element)->children.resize(3);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
			if (!Generate(node->getChild(1), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[1])) {
				return false;
			}
			if (!Generate(node->getChild(2), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[2])) {
				return false;
			}
		} else if (strcmp(node->getName(), "mm") == 0) {
			if (node->getNumChildren() != 4) {
				LOGERROR("Michaelis-menten function should have four parameters");
				return false;
			}
			*element = std::make_unique<SBMLRatelawElementFunctionMM>();
			(*element)->children.resize(4);
			if (!Generate(node->getChild(0), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[0])) {
				return false;
			}
			if (!Generate(node->getChild(1), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[1])) {
				return false;
			}
			if (!Generate(node->getChild(2), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[2])) {
				return false;
			}
			if (!Generate(node->getChild(3), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, forced_parameter_values, &(*element)->children[3])) {
				return false;
			}
		} else {
			LOGERROR("AST function with unknown name");
			return false;
		}
		break;

	default:
		LOGERROR("SBML AST node type %d not implemented", (int)node->getType());
		return false;
	}

	return true;
}

SBMLRatelawElement::~SBMLRatelawElement()
{
}

bool SBMLRatelawElement::ContainsSpeciesLookup(size_t species_ix)
{
	for (size_t i = 0; i < children.size(); i++) {
		if (children[i]->ContainsSpeciesLookup(species_ix)) {
			return true;
		}
	}
	return false;
}

Real SBMLRatelawElementPlus::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() > 1);
	Real result = children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
	for (size_t i = 1; i < children.size(); i++) {
		result += children[i]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
	}
	return result;
}

std::string SBMLRatelawElementPlus::GenerateEquation()
{
	ASSERT(children.size() > 1);
	std::string result = "(";
	result += children[0]->GenerateEquation();
	for (size_t i = 1; i < children.size(); i++) {
		result += "+";
		result += children[i]->GenerateEquation();
	}
	result += ")";
	return result;
}

std::string SBMLRatelawElementPlus::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() > 1);

	std::vector<std::string> children_result;
	bool nonzero = false;
	for (size_t i = 0; i < children.size(); i++) {
		std::string res = children[i]->GenerateDerivative(species_ix);
		children_result.push_back(res);
		if (!res.empty()) {
			nonzero = true;
		}
	}

	if (nonzero) {
		std::string result = "(";
		for (size_t i = 0; i < children.size(); i++) {
			if (!children_result[i].empty()) {
				result += "+";
				result += children_result[i];
			}
		}
		result += ")";
		return result;
	} else {
		return "";
	}
}

Real SBMLRatelawElementMinus::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 2);
	return children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters) - children[1]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
}

std::string SBMLRatelawElementMinus::GenerateEquation()
{
	ASSERT(children.size() == 2);
	std::string result = "(";
	result += children[0]->GenerateEquation();
	result += "-";
	result += children[1]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementMinus::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 2);

	std::string res1 = children[0]->GenerateDerivative(species_ix);
	std::string res2 = children[1]->GenerateDerivative(species_ix);

	if (res1.empty()) {
		if (res2.empty()) {
			return "";
		} else {
			std::string result = "(";
			result += "-";
			result += res2;
			result += ")";
			return result;
		}
	} else {
		if (res2.empty()) {
			return res1;
		} else {
			std::string result = "(";
			result += res1;
			result += "-";
			result += res2;
			result += ")";
			return result;
		}
	}
}

Real SBMLRatelawElementNegate::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 1);
	return -children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
}

std::string SBMLRatelawElementNegate::GenerateEquation()
{
	ASSERT(children.size() == 1);
	std::string result = "(";
	result += "-";
	result += children[0]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementNegate::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 1);
	std::string res = children[0]->GenerateDerivative(species_ix);
	if (res.empty()) {
		return "";
	} else {
		std::string result = "(";
		result += "-";
		result += res;
		result += ")";
		return result;
	}
}

Real SBMLRatelawElementTimes::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() > 1);
	Real result = children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
	for (size_t i = 1; i < children.size(); i++) {
		result *= children[i]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
	}
	return result;
}

std::string SBMLRatelawElementTimes::GenerateEquation()
{
	ASSERT(children.size() > 1);
	std::string result = "(";
	result += children[0]->GenerateEquation();
	for (size_t i = 1; i < children.size(); i++) {
		result += "*";
		result += children[i]->GenerateEquation();
	}
	result += ")";
	return result;
}

std::string SBMLRatelawElementTimes::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() > 1);

	std::vector<std::string> children_result;
	bool nonzero = false;
	for (size_t i = 0; i < children.size(); i++) {
		std::string res = children[i]->GenerateDerivative(species_ix);
		children_result.push_back(res);
		if (!res.empty()) {
			nonzero = true;
		}
	}

	if (nonzero) {
		// Product rule
		std::string result = "(";
		for (size_t i = 0; i < children.size(); i++) {
			if (!children_result[i].empty()) {
				for (size_t j = 0; j < children.size(); j++) {
					if (j != 0) {
						result += "*";
					}
					if (i == j) {
						result += children_result[i];
					} else {
						result += children[j]->GenerateEquation();
					}
				}
			}
		}
		result += ")";
		return result;
	} else {
		return "";
	}
}

Real SBMLRatelawElementDivide::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 2);
	return children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters) / children[1]->Evaluate(species, constant_species, parameters, non_sampled_parameters);
}

std::string SBMLRatelawElementDivide::GenerateEquation()
{
	ASSERT(children.size() == 2);
	std::string result = "(";
	result += children[0]->GenerateEquation();
	result += "/";
	result += children[1]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementDivide::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 2);

	std::string res1 = children[0]->GenerateDerivative(species_ix);
	std::string res2 = children[1]->GenerateDerivative(species_ix);
	std::string eq1 = children[0]->GenerateEquation();
	std::string eq2 = children[1]->GenerateEquation();

	if (res1.empty()) {
		if (res2.empty()) {
			return "";
		} else {
			std::string result = "(-";
			result += res2;
			result += "*";
			result += eq1;
			result += "/";
			result += "(";
			result += eq2;
			result += "*";
			result += eq2;
			result += ")";
			result += ")";
			return result;
		}
	} else {
		if (res2.empty()) {
			std::string result = "(";
			result += res1;
			result += "/";
			result += eq2;
			result += ")";
			return result;
		} else {
			std::string result = "(";
			result += "(";
			result += res1;
			result += "*";
			result += eq2;
			result += "-";
			result += res2;
			result += "*";
			result += eq1;
			result += ")";
			result += "/";
			result += "(";
			result += eq2;
			result += "*";
			result += eq2;
			result += ")";
			result += ")";
			return result;
		}
	}
}

Real SBMLRatelawElementLookupSpecies::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	return species[ix];
}

std::string SBMLRatelawElementLookupSpecies::GenerateEquation()
{
	std::string result = "species[";
	result += std::to_string((uint64)ix);
	result += "]";
	return result;
}

bool SBMLRatelawElementLookupSpecies::ContainsSpeciesLookup(size_t species_ix)
{
	return species_ix == ix;
}

std::string SBMLRatelawElementLookupSpecies::GenerateDerivative(size_t species_ix)
{
	if (species_ix == ix) {
#if ODE_SINGLE_PRECISION
		return "1.0f";
#else
		return "1.0";
#endif
	} else {
		return "";
	}
}

Real SBMLRatelawElementLookupConstantSpecies::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	return constant_species[ix];
}

std::string SBMLRatelawElementLookupConstantSpecies::GenerateEquation()
{
	std::string result = "constant_species[";
	result += std::to_string((uint64)ix);
	result += "]";
	return result;
}

bool SBMLRatelawElementLookupConstantSpecies::ContainsSpeciesLookup(size_t species_ix)
{
	return false;
}

std::string SBMLRatelawElementLookupConstantSpecies::GenerateDerivative(size_t species_ix)
{
	return "";
}

Real SBMLRatelawElementLookupParameter::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	return parameters[ix];
}

std::string SBMLRatelawElementLookupParameter::GenerateEquation()
{
	std::string result = "parameters[";
	result += std::to_string((uint64)ix);
	result += "]";
	return result;
}

bool SBMLRatelawElementLookupParameter::ContainsSpeciesLookup(size_t species_ix)
{
	return false;
}

std::string SBMLRatelawElementLookupParameter::GenerateDerivative(size_t species_ix)
{
	return "";
}

Real SBMLRatelawElementLookupNonSampledParameter::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	return non_sampled_parameters[ix];
}

std::string SBMLRatelawElementLookupNonSampledParameter::GenerateEquation()
{
	std::string result = "non_sampled_parameters[";
	result += std::to_string((uint64)ix);
	result += "]";
	return result;
}

bool SBMLRatelawElementLookupNonSampledParameter::ContainsSpeciesLookup(size_t species_ix)
{
	return false;
}

std::string SBMLRatelawElementLookupNonSampledParameter::GenerateDerivative(size_t species_ix)
{
	return "";
}

Real SBMLRatelawElementConstant::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	return constant;
}

std::string SBMLRatelawElementConstant::GenerateEquation()
{
#if ODE_SINGLE_PRECISION
	return std::to_string((long double)constant) + "f";
#else
	return std::to_string((long double)constant);
#endif
}

bool SBMLRatelawElementConstant::ContainsSpeciesLookup(size_t species_ix)
{
	return false;
}

std::string SBMLRatelawElementConstant::GenerateDerivative(size_t species_ix)
{
	return "";
}

Real SBMLRatelawElementFunctionExp::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 1);
	return exp(children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters));
}

std::string SBMLRatelawElementFunctionExp::GenerateEquation()
{
	ASSERT(children.size() == 1);
	std::string result = "exp(";
	result += children[0]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementFunctionExp::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 1);
	std::string res = children[0]->GenerateDerivative(species_ix);
	if (res.empty()) {
		return "";
	} else {
		std::string result = "exp(";
		result += children[0]->GenerateEquation();
		result += ")*(";
		result += res;
		result += ")";
		return result;
	}
}

Real SBMLRatelawElementFunctionLog::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 1);
	return log(children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters));
}

std::string SBMLRatelawElementFunctionLog::GenerateEquation()
{
	ASSERT(children.size() == 1);
	std::string result = "log(";
	result += children[0]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementFunctionLog::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 1);
	std::string res = children[0]->GenerateDerivative(species_ix);
	if (res.empty()) {
		return "";
	} else {
		std::string result = "(double)1.0/(";
		result += children[0]->GenerateEquation();
		result += ")";
		return result;
	}
}

Real SBMLRatelawElementFunctionPow::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 2);
	return safepow(children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters), children[1]->Evaluate(species, constant_species, parameters, non_sampled_parameters));
}

std::string SBMLRatelawElementFunctionPow::GenerateEquation()
{
	std::string result = "safepow(";
	result += children[0]->GenerateEquation();
	result += ",";
	result += children[1]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementFunctionPow::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 2);

	std::string res1 = children[0]->GenerateDerivative(species_ix);
	std::string res2 = children[1]->GenerateDerivative(species_ix);
	if (!res2.empty()) {
		LOGERROR("Not implemented");
	}

	std::string result;
	if (!res1.empty()) {
		result = "(";
		result += res1;
		result += "*pow(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
#if ODE_SINGLE_PRECISION
		result += "-1.0)*(";
#else
		result += "-1.0f)*(";
#endif
		result += children[1]->GenerateEquation();
		result += "))";
	}
	return result;
}

Real SBMLRatelawElementFunctionHill::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 3);
	return hill_function(children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters),
						 children[1]->Evaluate(species, constant_species, parameters, non_sampled_parameters),
						 children[2]->Evaluate(species, constant_species, parameters, non_sampled_parameters));
}

std::string SBMLRatelawElementFunctionHill::GenerateEquation()
{
	std::string result;
	std::string n = children[2]->GenerateEquation();
#if USE_INTEGER_HILL_OPTIMIZATIONS
	if (n == "1.000000" || n == "1.000000f") {
		result = "hill_function_fixedn1(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ")";
	} else if (n == "2.000000" || n == "2.000000f") {
		result = "hill_function_fixedn2(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ")";
	} else if (n == "4.000000" || n == "4.000000f") {
		result = "hill_function_fixedn4(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ")";
	} else if (n == "10.000000" || n == "10.000000f") {
		result = "hill_function_fixedn10(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ")";
	} else if (n == "16.000000" || n == "16.000000f") {
		result = "hill_function_fixedn16(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ")";
	} else if (n == "100.000000" || n == "100.000000f") {
		result = "hill_function_fixedn100(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ")";
	} else {
#endif
		result = "hill_function(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ",";
		result += children[2]->GenerateEquation();
		result += ")";
#if USE_INTEGER_HILL_OPTIMIZATIONS
	}
#endif
	return result;
}

std::string SBMLRatelawElementFunctionHill::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 3);
	std::string result;
	SBMLRatelawElementLookupSpecies* lookup_species = dynamic_cast<SBMLRatelawElementLookupSpecies*>(children[0].get());
	if (lookup_species && lookup_species->ix == species_ix) {
		std::string n = children[2]->GenerateEquation();
#if USE_INTEGER_HILL_OPTIMIZATIONS
		if (n == "1.000000" || n == "1.000000f") {
			result = "hill_function_derivative_fixedn1(";
			result += children[0]->GenerateEquation();
			result += ",";
			result += children[1]->GenerateEquation();
			result += ")";
		} else {
#endif
			result = "hill_function_derivative(";
			result += children[0]->GenerateEquation();
			result += ",";
			result += children[1]->GenerateEquation();
			result += ",";
			result += n;
			result += ")";
#if USE_INTEGER_HILL_OPTIMIZATIONS
		}
#endif
	} else {
		for (size_t i = 0; i < children.size(); i++) {
			std::string res = children[i]->GenerateDerivative(species_ix);
			if (!res.empty()) {
				if (i == 0) {
					std::string n = children[2]->GenerateEquation();
#if USE_INTEGER_HILL_OPTIMIZATIONS
					if (n == "1.000000" || n == "1.000000f") {
						result = "hill_function_derivative_fixedn1(";
						result += children[0]->GenerateEquation();
						result += ",";
						result += children[1]->GenerateEquation();
						result += ")*(";
						result += res;
						result += ")";
					} else {
#endif
						result = "hill_function_derivative(";
						result += children[0]->GenerateEquation();
						result += ",";
						result += children[1]->GenerateEquation();
						result += ",";
						result += n;
						result += ")*(";
						result += res;
						result += ")";
#if USE_INTEGER_HILL_OPTIMIZATIONS
					}
#endif
				} else {
					LOGERROR("Not implemented");
				}
			}
		}
	}
	return result;
}

Real SBMLRatelawElementFunctionMM::Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)
{
	ASSERT(children.size() == 4);
	return michaelis_menten_function(
		children[0]->Evaluate(species, constant_species, parameters, non_sampled_parameters),
		children[1]->Evaluate(species, constant_species, parameters, non_sampled_parameters),
		children[2]->Evaluate(species, constant_species, parameters, non_sampled_parameters),
		children[3]->Evaluate(species, constant_species, parameters, non_sampled_parameters)
	);
}

std::string SBMLRatelawElementFunctionMM::GenerateEquation()
{
	std::string result;
	result = "michaelis_menten_function(";
	result += children[0]->GenerateEquation();
	result += ",";
	result += children[1]->GenerateEquation();
	result += ",";
	result += children[2]->GenerateEquation();
	result += ",";
	result += children[3]->GenerateEquation();
	result += ")";
	return result;
}

std::string SBMLRatelawElementFunctionMM::GenerateDerivative(size_t species_ix)
{
	ASSERT(children.size() == 4);
	std::string result;

	for (size_t i = 0; i < 2; i++) {
		std::string child_deriv = children[i]->GenerateDerivative(species_ix);
		if (!child_deriv.empty()) {
			LOGERROR("Not implemented");
		}
	}

	std::string enzyme_deriv = children[2]->GenerateDerivative(species_ix);
	if (!enzyme_deriv.empty()) {
		result = "michaelis_menten_derivative_enzyme(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ",";
		result += children[2]->GenerateEquation();
		result += ",";
		result += children[3]->GenerateEquation();
		result += ")*(";
		result += enzyme_deriv;
		result += ")";
	}

	std::string substrate_deriv = children[3]->GenerateDerivative(species_ix);
	if (!substrate_deriv.empty()) {
		result = "michaelis_menten_derivative_substrate(";
		result += children[0]->GenerateEquation();
		result += ",";
		result += children[1]->GenerateEquation();
		result += ",";
		result += children[2]->GenerateEquation();
		result += ",";
		result += children[3]->GenerateEquation();
		result += ")*(";
		result += substrate_deriv;
		result += ")";
	}

	if (!enzyme_deriv.empty() && !substrate_deriv.empty()) {
		LOGERROR("Not implemented");
	}

	return result;
}
