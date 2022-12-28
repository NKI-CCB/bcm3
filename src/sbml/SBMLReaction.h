#pragma once

#include <sbml/SBMLTypes.h>
#include "VariableSet.h"

class SBMLModel;
class SBMLSpecies;
class SBMLRatelawElement;

class SBMLReaction
{
public:
	SBMLReaction();
	~SBMLReaction();

	enum EType {
		Association,
		Dissociation,
		Transition,
		CatalyzedTransition,
		Transcription,
		Translation,
		PhysicalStimulation,

		Invalid
	};

	bool Initialize(const Reaction* reaction, const SBMLModel* model);
	bool SetSpeciesVector(const std::vector<std::string>& species_ids, const std::vector<std::string>& constant_species_ids);
	bool SetVariableSet(const bcm3::VariableSet* variables);
	bool SetNonSampledParameters(const std::vector<std::string>& parameter_names);
	bool AddRates(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters, OdeReal* rates);
	std::string GetRateEquation(size_t species_ix);
	void AddRateDerivativeEquation(size_t species_ix, size_t deriv_ix, std::string& eqn);
	bool PostInitialize(const Model* model, const std::map<std::string, Real>& fixed_parameter_values);

	inline const std::string& GetId() const { return id; }
	inline const std::vector<std::string>& GetReactants() const { return reactants; }
	inline const std::vector<std::string>& GetProducts() const { return products; }
	inline Real GetReactantStoichiometry(size_t ix) const { return reactant_stoichiometry[ix]; }
	inline Real GetProductStoichiometry(size_t ix) const { return product_stoichiometry[ix]; }

private:
	void MapRateLawSpecies(const std::vector<std::string>& species_ids, const ASTNode* node);
	void MapRateLawConstantSpecies(const std::vector<std::string>& constant_species_ids, const ASTNode* node);
	void MapRateLawParameters(const bcm3::VariableSet* variables, const ASTNode* node);
	void MapRateLawNonSampledParameters(const std::vector<std::string>& parameter_names, const ASTNode* node);

	std::string id;
	EType type;
	std::vector<std::string> reactants;
	std::vector<std::string> products;
	std::vector<size_t> reactant_indices;
	std::vector<size_t> product_indices;
	std::vector<Real> reactant_stoichiometry;
	std::vector<Real> product_stoichiometry;
	std::map<std::string, size_t> species_index_map;
	std::map<std::string, size_t> constant_species_index_map;
	std::map<std::string, size_t> parameter_index_map;
	std::map<std::string, size_t> non_sampled_parameter_index_map;

	std::unique_ptr<ASTNode> rate_law;
	std::unique_ptr<SBMLRatelawElement> evaluate_rate_law;
};
