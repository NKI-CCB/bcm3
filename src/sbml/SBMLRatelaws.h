#pragma once

#include <sbml/SBMLTypes.h>
#include "LinearAlgebraSelector.h"

class SBMLModel;
class SBMLSpecies;

class SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElement();
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters) = 0;
	virtual std::string GenerateEquation() = 0;
	virtual bool ContainsSpeciesLookup(size_t species_ix);
	virtual std::string GenerateDerivative(size_t species_ix) = 0;
	std::vector< std::unique_ptr<SBMLRatelawElement> > children;

	static bool Generate(const ASTNode* node, const Model* model, const std::map<std::string, size_t>& species_index_map,
																  const std::map<std::string, size_t>& parameter_index_map,
															      const std::map<std::string, size_t>& constant_species_index_map,
																  const std::map<std::string, size_t>& non_sampled_parameter_index_map,
																  const std::map<std::string, Real>& forced_parameter_values,
																  std::unique_ptr<SBMLRatelawElement>* element);
};

class SBMLRatelawElementPlus : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementPlus() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementMinus : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementMinus() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementNegate : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementNegate() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementTimes : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementTimes() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementDivide : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementDivide() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementLookupSpecies : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementLookupSpecies() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual bool ContainsSpeciesLookup(size_t species_ix);
	virtual std::string GenerateDerivative(size_t species_ix);
	size_t ix;
};

class SBMLRatelawElementLookupConstantSpecies : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementLookupConstantSpecies() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual bool ContainsSpeciesLookup(size_t species_ix);
	virtual std::string GenerateDerivative(size_t species_ix);
	size_t ix;
};

class SBMLRatelawElementLookupParameter : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementLookupParameter() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual bool ContainsSpeciesLookup(size_t species_ix);
	virtual std::string GenerateDerivative(size_t species_ix);
	size_t ix;
};

class SBMLRatelawElementLookupNonSampledParameter : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementLookupNonSampledParameter() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual bool ContainsSpeciesLookup(size_t species_ix);
	virtual std::string GenerateDerivative(size_t species_ix);
	size_t ix;
};

class SBMLRatelawElementConstant : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementConstant() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual bool ContainsSpeciesLookup(size_t species_ix);
	virtual std::string GenerateDerivative(size_t species_ix);
	Real constant;
};

class SBMLRatelawElementFunctionExp : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementFunctionExp() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementFunctionLog : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementFunctionLog() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementFunctionPow : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementFunctionPow() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementFunctionHill : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementFunctionHill() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};

class SBMLRatelawElementFunctionMM : public SBMLRatelawElement
{
public:
	virtual ~SBMLRatelawElementFunctionMM() {}
	virtual Real Evaluate(const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	virtual std::string GenerateEquation();
	virtual std::string GenerateDerivative(size_t species_ix);
};
