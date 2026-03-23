#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "VariabilityDescription.h"

VariabilityDescription::VariabilityDescription()
	: distribution(EDistribution::Invalid)
{

}

VariabilityDescription::~VariabilityDescription()
{

}

std::unique_ptr<VariabilityDescription> VariabilityDescription::Create(const boost::property_tree::ptree& xml_node)
{
	std::unique_ptr<VariabilityDescription> desc;

	try {
		desc = std::make_unique<VariabilityDescription>();
		if (!desc->Load(xml_node)) {
			desc.reset();
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		desc.reset();
	}

	return desc;
}

bool VariabilityDescription::PostInitialize(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model)
{
    for (const auto& it : variables) {
        if (!it->PostInitialize(varset, non_sampled_parameter_names, model)) {
            return false;
        }
    }

    return true;
}

VectorReal VariabilityDescription::GetPseudorandomVector(const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	size_t D = variables.size();

	switch (distribution) {
	case EDistribution::DiagonalGaussian:
		{
			VectorReal v(D);
			for (size_t i = 0; i < D; ++i) {
				v(i) = bcm3::QuantileNormal(sobol_sequence(sobol_sequence_ix++), 0, 1);
				v(i) *= exp(value_references[i].GetValue(transformed_values, non_sampled_parameters));
			}
			return v;
		}

	case EDistribution::FullGaussian:
		{
			// Log-cholesky parametrization as described by Pinheiro & Bates, Unconstrained parametrizations for variance-covariance matrices, Statistics and Computing 1996
			// The diagonal terms are on log scale, the off-diagonal on natural scale
			MatrixReal cholesky_L(D, D);
			int value_reference_i = 0;
			for (int i = 0; i < D; i++) {
				cholesky_L(i, i) = exp(value_references[value_reference_i++].GetValue(transformed_values, non_sampled_parameters));
				for (int j = 0; j < i; j++) {
					cholesky_L(i, j) = value_references[value_reference_i++].GetValue(transformed_values, non_sampled_parameters);
				}
				for (int j = i + 1; j < D; j++) {
					cholesky_L(i, j) = 0.0;
				}
			}

			// Normally distributed pseudorandom values
			VectorReal v(D);
			for (int i = 0; i < D; i++) {
				v(i) = bcm3::QuantileNormal(sobol_sequence(sobol_sequence_ix++), 0, 1);
			}

			// Transform by covariance matrix
			v = cholesky_L * v;

			return v;
		}

	default:
		ASSERT(false);
		break;
	}

	return VectorReal();
}

void VariabilityDescription::ApplyVariabilityEntryTime(Real& value, const VectorReal& pseudorandom_vector, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	ASSERT(pseudorandom_vector.size() == variables.size());
	for (int i = 0; i < variables.size(); i++) {
		variables[i]->ApplyVariabilityEntryTime(value, pseudorandom_vector(i), transformed_values, non_sampled_parameters, is_initial_cell);
	}
}

void VariabilityDescription::ApplyVariabilityParameter(const std::string& parameter, OdeReal& value, const VectorReal& pseudorandom_vector, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	ASSERT(pseudorandom_vector.size() == variables.size());
	for (int i = 0; i < variables.size(); i++) {
		variables[i]->ApplyVariabilityEntryTime(value, pseudorandom_vector(i), transformed_values, non_sampled_parameters, is_initial_cell);
	}
}

void VariabilityDescription::ApplyVariabilitySpecies(const std::string& species, OdeReal& value, const VectorReal& pseudorandom_vector, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	ASSERT(pseudorandom_vector.size() == variables.size());
	for (int i = 0; i < variables.size(); i++) {
		variables[i]->ApplyVariabilityEntryTime(value, pseudorandom_vector(i), transformed_values, non_sampled_parameters, is_initial_cell);
	}
}

size_t VariabilityDescription::GetNumDimensions() const
{
	return variables.size();
}

bool VariabilityDescription::Load(const boost::property_tree::ptree& xml_node)
{
	try {
		BOOST_FOREACH(const boost::property_tree::ptree::value_type & var, xml_node.get_child("")) {
			if (var.first == "variable") {
				std::unique_ptr< VariabilityDescriptionVariable> vdv = VariabilityDescriptionVariable::Create(var.second);
				if (vdv) {
					variables.push_back(std::move(vdv));
				} else {
					return false;
				}
			}
		}
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	std::string distribution_str = xml_node.get<std::string>("<xmlattr>.distribution");
	if (distribution_str == "full_gaussian") {
		distribution = EDistribution::FullGaussian;
	} else if (distribution_str == "diagonal_gaussian") {
		distribution = EDistribution::DiagonalGaussian;
	} else {
		LOGERROR("Unknown distribution \"%s\" in variability description", distribution_str.c_str());
		return false;
	}

	switch (distribution) {
	case EDistribution::DiagonalGaussian:
		{
			for (size_t i = 0; i < variables.size(); ++i) {
				value_references.push_back(variables[i]->scale);
			}
		}
		break;

	case EDistribution::FullGaussian:
		{
		for (size_t i = 0; i < variables.size(); ++i) {
			value_references.push_back(variables[i]->scale);
			for (int j = 0; j < i; j++) {

			}
		}
		}
		break;
	}

	return true;
}
