#include "Utils.h"
#include "SBMLModel.h"
#include "SBMLRatelaws.h"
#include "SBMLReaction.h"
#include "SBMLSpecies.h"

SBMLReaction::SBMLReaction()
	: type(Invalid)
{
}

SBMLReaction::~SBMLReaction()
{
}

bool SBMLReaction::Initialize(const Reaction* reaction, const SBMLModel* model)
{
	id = reaction->getId();

	if (reaction->getKineticLaw() == NULL) {
		LOGERROR("Reaction \"%s\" does not have a kinetic law", id.c_str());
		return false;
	}

	rate_law = std::move(std::unique_ptr<ASTNode>(reaction->getKineticLaw()->getMath()->deepCopy()));

	const XMLNode* annotation = reaction->getAnnotation();
	if (annotation) {
		for (unsigned int i = 0; i < annotation->getNumChildren(); i++) {
			const XMLNode& child = annotation->getChild(i);
			if (child.getPrefix() == "celldesigner" && child.getName() == "extension") {
				std::string reaction_type_str = child.getChild("reactionType").getChild(0).getCharacters();
				if (reaction_type_str.empty()) {
					LOGERROR("Can't find reaction type for reaction %s - is the SBML file a CellDesigner file?", id.c_str());
					return false;
				}
			
				if (reaction_type_str == "HETERODIMER_ASSOCIATION") {
					type = Association;
				} else if (reaction_type_str == "DISSOCIATION") {
					type = Dissociation;
				} else if (reaction_type_str == "STATE_TRANSITION") {
					type = Transition;
				} else if (reaction_type_str == "TRANSCRIPTION") {
					type = Transcription;
				} else if (reaction_type_str == "TRANSLATION") {
					type = Translation;
				} else if (reaction_type_str == "PHYSICAL_STIMULATION") {
					type = PhysicalStimulation;
				} else {
					LOGERROR("Unrecognized reaction type %s for reaction %s", reaction_type_str.c_str(), id.c_str());
					return false;
				}
			}
		}
	}
	
	for (unsigned int i = 0; i < reaction->getNumReactants(); i++) {
		const std::string& species_id = reaction->getReactant(i)->getSpecies();
		const SBMLSpecies* species = model->GetSpecies(species_id);
		if (species && species->GetType() != SBMLSpecies::Sink) {
			reactants.push_back(species_id);
			reactant_stoichiometry.push_back(reaction->getReactant(i)->getStoichiometry());
		}
	}

	for (unsigned int i = 0; i < reaction->getNumProducts(); i++) {
		const std::string& species_id = reaction->getProduct(i)->getSpecies();
		const SBMLSpecies* species = model->GetSpecies(species_id);
		if (species && species->GetType() != SBMLSpecies::Sink) {
			products.push_back(species_id);
			product_stoichiometry.push_back(reaction->getProduct(i)->getStoichiometry());
		}
	}

	return true;
}

bool SBMLReaction::SetSpeciesVector(const std::vector<std::string>& species_ids, const std::vector<std::string>& constant_species_ids)
{
	// Check if any of the named AST nodes map to a species
	MapRateLawSpecies(species_ids, rate_law.get());
	MapRateLawConstantSpecies(constant_species_ids, rate_law.get());

	// Map reactants and products
	for (std::vector<std::string>::iterator ri = reactants.begin(); ri != reactants.end(); ++ri) {
		bool found = false;
		for (size_t i = 0; i < species_ids.size(); i++) {
			if (species_ids[i] == *ri) {
				reactant_indices.push_back(i);
				found = true;
				break;
			}
		}

		if (!found) {
			LOGERROR("Can't find reactant \"%s\" in the species vector for reaction %s", ri->c_str(), id.c_str());
			return false;
		}
	}
	
	for (std::vector<std::string>::iterator pi = products.begin(); pi != products.end(); ++pi) {
		bool found = false;
		for (size_t i = 0; i < species_ids.size(); i++) {
			if (species_ids[i] == *pi) {
				product_indices.push_back(i);
				found = true;
				break;
			}
		}

		if (!found) {
			LOGERROR("Can't find product %s in the species vector", pi->c_str());
			return false;
		}
	}

	return true;
}

bool SBMLReaction::SetVariableSet(const bcm3::VariableSet* variables)
{
	// Check if any of the named AST nodes map to a variable
	MapRateLawParameters(variables, rate_law.get());
	return true;
}

bool SBMLReaction::SetNonSampledParameters(const std::vector<std::string>& parameter_names)
{
	// Check if any of the named AST nodes map to a variable
	MapRateLawNonSampledParameters(parameter_names, rate_law.get());
	return true;
}

bool SBMLReaction::PostInitialize(const Model* model, const std::map<std::string, Real>& fixed_parameter_values)
{
	evaluate_rate_law.reset();
	return SBMLRatelawElement::Generate(rate_law.get(), model, species_index_map, parameter_index_map, constant_species_index_map, non_sampled_parameter_index_map, fixed_parameter_values, &evaluate_rate_law);
}

bool SBMLReaction::AddRates(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters, OdeReal* rates)
{
	Real r = evaluate_rate_law->Evaluate(species, constant_species, parameters, non_sampled_parameters);
	for (std::vector<size_t>::iterator ri = reactant_indices.begin(); ri != reactant_indices.end(); ++ri) {
		rates[*ri] -= r;
	}
	for (std::vector<size_t>::iterator pi = product_indices.begin(); pi != product_indices.end(); ++pi) {
		rates[*pi] += r;
	}
	return true;
}

std::string SBMLReaction::GetRateEquation(size_t species_ix)
{
	return evaluate_rate_law->GenerateEquation();
}

void SBMLReaction::AddRateDerivativeEquation(size_t species_ix, size_t deriv_ix, std::string& eqn)
{
	for (size_t ri = 0; ri < reactant_indices.size(); ri++) {
		if (reactant_indices[ri] == species_ix) {
			std::string derivative = evaluate_rate_law->GenerateDerivative(deriv_ix);
			if (!derivative.empty()) {
				eqn += "-";
				if (reactant_stoichiometry[ri] != 0.0) {
					if (reactant_stoichiometry[ri] != 1.0) {
#if ODE_SINGLE_PRECISION
						eqn += std::to_string(reactant_stoichiometry[ri]) + "f*";
#else
						eqn += std::to_string(reactant_stoichiometry[ri]) + "*";
#endif
					}
					eqn += derivative;
				}
			}
		}
	}
	for (size_t pi = 0; pi < product_indices.size(); pi++) {
		if (product_indices[pi] == species_ix) {
			std::string derivative = evaluate_rate_law->GenerateDerivative(deriv_ix);
			if (!derivative.empty()) {
				eqn += "+";
				if (product_stoichiometry[pi] != 0.0) {
					if (product_stoichiometry[pi] != 1.0) {
#if ODE_SINGLE_PRECISION
						eqn += std::to_string(product_stoichiometry[pi]) + "f*";
#else
						eqn += std::to_string(product_stoichiometry[pi]) + "*";
#endif
					}
					eqn += derivative;
				}
			}
		}
	}
}

void SBMLReaction::MapRateLawSpecies(const std::vector<std::string>& species_ids, const ASTNode* node)
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

void SBMLReaction::MapRateLawConstantSpecies(const std::vector<std::string>& constant_species_ids, const ASTNode* node)
{
	if (node->getType() == AST_NAME) {
		const char* name = node->getName();
		for (size_t i = 0; i < constant_species_ids.size(); i++) {
			if (constant_species_ids[i] == name) {
				constant_species_index_map[name] = i;
			}
		}
	} else {
		for (unsigned int i = 0; i < node->getNumChildren(); i++) {
			MapRateLawConstantSpecies(constant_species_ids, node->getChild(i));
		}
	}
}

void SBMLReaction::MapRateLawParameters(const bcm3::VariableSet* variables, const ASTNode* node)
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

void SBMLReaction::MapRateLawNonSampledParameters(const std::vector<std::string>& parameter_names, const ASTNode* node)
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