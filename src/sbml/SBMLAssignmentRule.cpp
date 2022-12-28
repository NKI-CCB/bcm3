#include "Utils.h"
#include "SBMLAssignmentRule.h"

SBMLAssignmentRule::SBMLAssignmentRule()
{
}

SBMLAssignmentRule::~SBMLAssignmentRule()
{
}
	
bool SBMLAssignmentRule::Initialize(const Rule* rule, const SBMLModel* model)
{
	if (!rule->isAssignment()) {
		LOGERROR("SBML rule is not an assignment");
		return false;
	}

	rate_law = std::move(std::unique_ptr<ASTNode>(rule->getMath()->deepCopy()));
	TargetVariable = rule->getVariable();

	return true;
}

bool SBMLAssignmentRule::SetSpeciesVector(const std::vector<std::string>& species_ids)
{
	// Check if any of the named AST nodes map to a species
	MapRateLawSpecies(species_ids, rate_law.get());
	return true;
}

bool SBMLAssignmentRule::SetVariableSet(const bcm3::VariableSet* variables)
{
	// Check if any of the named AST nodes map to a variable
	MapRateLawParameters(variables, rate_law.get());
	return true;
}

bool SBMLAssignmentRule::SetNonSampledParameters(const std::vector<std::string>& parameter_names)
{
	// Check if any of the named AST nodes map to a variable
	MapRateLawNonSampledParameters(parameter_names, rate_law.get());
	return true;
}

bool SBMLAssignmentRule::PostInitialize(const Model* model, const std::map<std::string, Real>& fixed_parameter_values)
{
	evaluate_rate_law.reset();
	return SBMLRatelawElement::Generate(rate_law.get(), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, fixed_parameter_values, &evaluate_rate_law);
}

bool SBMLAssignmentRule::Calculate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters, OdeReal* value) const
{
	ASSERT(value);
	*value = evaluate_rate_law->Evaluate(species, constant_species, parameters, non_sampled_parameters);
	return true;
}

void SBMLAssignmentRule::MapRateLawSpecies(const std::vector<std::string>& species_ids, const ASTNode* node)
{
	if (node->getType() == AST_NAME) {
		const char* name = node->getName();
		for (size_t i = 0; i < species_ids.size(); i++) {
			if (species_ids[i] == name) {
				species_index_map[name] = i;
			}
		}
	} else {
		for (unsigned int i = 0; i < node->getNumChildren(); i++) {
			MapRateLawSpecies(species_ids, node->getChild(i));
		}
	}
}

void SBMLAssignmentRule::MapRateLawParameters(const bcm3::VariableSet* variables, const ASTNode* node)
{
	if (node->getType() == AST_NAME) {
		const char* name = node->getName();
		for (size_t i = 0; i < variables->GetNumVariables(); i++) {
			if (variables->GetVariableName(i) == name) {
				parameter_index_map[name] = i;
			}
		}
	} else {
		for (unsigned int i = 0; i < node->getNumChildren(); i++) {
			MapRateLawParameters(variables, node->getChild(i));
		}
	}
}

void SBMLAssignmentRule::MapRateLawNonSampledParameters(const std::vector<std::string>& parameter_names, const ASTNode* node)
{
	if (node->getType() == AST_NAME) {
		const char* name = node->getName();
		for (size_t i = 0; i < parameter_names.size(); i++) {
			if (parameter_names[i] == name) {
				non_sampled_parameter_index_map[name] = i;
			}
		}
	} else {
		for (unsigned int i = 0; i < node->getNumChildren(); i++) {
			MapRateLawNonSampledParameters(parameter_names, node->getChild(i));
		}
	}
}