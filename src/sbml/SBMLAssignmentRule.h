#pragma once

#include <sbml/SBMLTypes.h>
#include "SBMLRatelaws.h"
#include "VariableSet.h"

class SBMLAssignmentRule
{
public:
	SBMLAssignmentRule();
	~SBMLAssignmentRule();
	
	bool Initialize(const Rule* rule, const SBMLModel* model);
	bool SetSpeciesVector(const std::vector<std::string>& species_ids);
	bool SetVariableSet(const bcm3::VariableSet* variables);
	bool SetNonSampledParameters(const std::vector<std::string>& parameter_names);
	bool PostInitialize(const Model* model, const std::map<std::string, Real>& fixed_parameter_values);

	bool Calculate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters, OdeReal* value) const;
	inline const std::string& GetVariable() const { return TargetVariable; }

private:
	void MapRateLawSpecies(const std::vector<std::string>& species_ids, const ASTNode* node);
	void MapRateLawParameters(const bcm3::VariableSet* variables, const ASTNode* node);
	void MapRateLawNonSampledParameters(const std::vector<std::string>& parameter_names, const ASTNode* node);

	std::string TargetVariable;
	std::map<std::string, size_t> species_index_map;
	std::map<std::string, size_t> constant_species_index_map;
	std::map<std::string, size_t> parameter_index_map;
	std::map<std::string, size_t> non_sampled_parameter_index_map;
	std::unique_ptr<ASTNode> rate_law;
	std::unique_ptr<SBMLRatelawElement> evaluate_rate_law;
};
