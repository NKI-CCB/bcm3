#pragma once

#define BCM3_SBML_INCLUDE_SOLVERS false

class ODESolverCVODE;
class ExperimentalCondition;
class SBMLAssignmentRule;
class SBMLReaction;
class SBMLSpecies;

#include <sbml/SBMLTypes.h>
#include "VariableSet.h"
#include "ODESolverCVODE.h"
#include "SBMLModelParameters.h"

class SBMLModel
{
public:
	SBMLModel(size_t numthreads);
	~SBMLModel();

	bool LoadSBML(const std::string& fn);
	bool SetVariableSet(std::shared_ptr<const bcm3::VariableSet> variables);
	bool AddNonSampledParameters(const std::vector<std::string>& non_sampled_parameters);
	bool AddFixedParameter(const std::string& parameter, Real value);
	bool Initialize();

	std::string GenerateCode();
	std::string GenerateJacobianCode(size_t i, size_t j);
	bool CalculateDerivativePublic(OdeReal t, const OdeReal* y, OdeReal* dydt, const OdeReal* constant_species_y, const OdeReal* transformed_variables, const OdeReal* non_sampled_parameters) const;
	
#if BCM3_SBML_INCLUDE_SOLVERS
	void DumpODEStatistics(const std::string& base_filename);
	void ResetForcedVariables(size_t threadix);
	bool SetForcedVariable(const std::string& variable_name, Real value, size_t threadix);
	bool Simulate(const Real* variables, const Real* initial_conditions, const VectorReal& timepoints, size_t threadix);
	bool GetOutput(size_t species_ix, VectorReal& output, size_t threadix);
	bool GetOutput(const std::string& species_full_name, VectorReal& output, size_t threadix);
#endif

	void GetParameters(std::vector<std::string>& parameters) const;
	size_t GetNumSimulatedSpecies() const;
	size_t GetNumCVodeSpecies() const;
	size_t GetNumConstantSpecies() const;
	const SBMLSpecies* GetSimulatedSpecies(size_t i) const;
	const SBMLSpecies* GetCVodeSpecies(size_t i) const;
	const SBMLSpecies* GetConstantSpecies(size_t i) const;
	const std::string& GetSimulatedSpeciesName(size_t i) const;
	const std::string& GetCVodeSpeciesName(size_t i) const;
	const std::string& GetConstantSpeciesName(size_t i) const;
	std::string GetSimulatedSpeciesFullName(size_t i) const;
	size_t GetSimulatedSpeciesByName(const std::string& name, bool log_error = true) const;
	inline size_t GetSimulatedSpeciesFromCVodeSpecies(size_t cvode_i) const { return cvode_to_simulated_map[cvode_i]; }
	inline size_t GetSimulatedSpeciesFromConstantSpecies(size_t constant_i) const { return constant_to_simulated_map[constant_i]; }
	size_t GetCVodeSpeciesByName(const std::string& name, bool log_error = true) const;
	size_t GetConstantSpeciesByName(const std::string& name, bool log_error = true) const;
	void GetSimulatedSpeciesFullNames(std::vector<std::string>& species_names) const;
	const SBMLSpecies* GetSpecies(const std::string& id) const;
	const SBMLSpecies* GetSpeciesFromFullName(const std::string& id) const;
	inline size_t GetNumAssignmentRules() const { return AssignmentRules.size(); }
	inline const SBMLAssignmentRule& GetAssignmentRule(size_t i) const { return *AssignmentRules[i]; }

private:
#if BCM3_SBML_INCLUDE_SOLVERS
	bool CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user) const;
	bool CalculateJacobian(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeReal** jac, void* user) const;
	void CalculateAssignments(const OdeReal* y, const OdeReal* variables, const OdeReal* constant_species_y, const OdeReal* non_sampled_parameters, OdeReal* out) const;
#endif

	SBMLDocument* document;
	const Model* model;
	std::shared_ptr<const bcm3::VariableSet> variables;
	
	std::map< std::string, std::unique_ptr<SBMLSpecies> > Species;
	std::map< std::string, std::unique_ptr<SBMLReaction> > Reactions;
	std::vector<std::string> SimulatedSpecies;
	std::vector<std::string> CVodeSpecies;
	std::vector<std::string> ConstantSpecies;
	std::vector< std::unique_ptr<SBMLAssignmentRule> > AssignmentRules;
	std::vector<size_t> AssignmentRulesTargetIx;
	std::vector<bool> SpeciesIsAssigned;
	std::vector<size_t> cvode_to_simulated_map;
	std::vector<size_t> constant_to_simulated_map;
	std::map<std::string, Real> fixed_parameter_values;

#if BCM3_SBML_INCLUDE_SOLVERS
	struct ParallelSolver {
		ParallelSolver() : variables(NULL) {}
		std::unique_ptr<ODESolverCVODE> solver;
		std::unique_ptr<SBMLModelParameters> parameters;
		std::vector<OdeReal> y;
		std::vector<OdeReal> constant_species_y;
		const Real* variables;
		VectorReal transformed_variables;
		OdeVectorReal non_sampled_variables;
		OdeVectorReal transformed_variables_use_in_ode;
		OdeMatrixReal cvode_trajectories;
		MatrixReal trajectories;
	};
	std::vector<ParallelSolver> Solvers;
#endif
};
