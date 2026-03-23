#pragma once

#include "SBMLModel.h"
#include "VariableSet.h"
#include "VariabilityDescriptionVariable.h"

#include <boost/random/sobol.hpp>

class VariabilityDescription
{
public:
	VariabilityDescription();
	~VariabilityDescription();

	static std::unique_ptr<VariabilityDescription> Create(const boost::property_tree::ptree& xml_node);
	bool PostInitialize(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model);

	VectorReal GetPseudorandomVector(const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const;

	void ApplyVariabilityEntryTime(Real& value, const VectorReal& pseudorandom_vector, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;
	void ApplyVariabilityParameter(const std::string& parameter, OdeReal& value, const VectorReal& pseudorandom_vector, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;
	void ApplyVariabilityInitialCondition(const std::string& species, OdeReal& value, const VectorReal& pseudorandom_vector, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;

	size_t GetNumDimensions() const;

private:
	bool Load(const boost::property_tree::ptree& xml_node);

	enum class EDistribution {
		FullGaussian,
		DiagonalGaussian,
		Invalid
	};

	std::vector< std::unique_ptr<VariabilityDescriptionVariable> > variables;
	EDistribution distribution;
	std::vector<ValueReference> covariance_values;
};
