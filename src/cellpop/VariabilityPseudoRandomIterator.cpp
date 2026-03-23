#include "Utils.h"
#include "VariabilityPseudoRandomIterator.h"

#include <boost/random/uniform_01.hpp>

VariabilityPseudoRandomIterator::VariabilityPseudoRandomIterator()
{ 
}

VariabilityPseudoRandomIterator::~VariabilityPseudoRandomIterator()
{
}

void VariabilityPseudoRandomIterator::Initialize(size_t dimensions, size_t initial_number_of_cells, size_t max_number_of_cells)
{
	sobol_sequence = std::make_shared< boost::random::sobol >(dimensions);
	sobol_sequence_values.resize(initial_number_of_cells * 100);
	boost::random::uniform_01<Real> unif;
	for (int i = 0; i < sobol_sequence_values.size(); i++) {
		sobol_sequence_values[i].resize(dimensions);
		for (int j = 0; j < dimensions; j++) {
			sobol_sequence_values[i](j) = unif(*sobol_sequence);
		}
	}
	sobol_sequence_indices.resize(max_number_of_cells, std::numeric_limits<size_t>::max());
}

void VariabilityPseudoRandomIterator::SetIndex(size_t cell_ix, size_t variability_index)
{
	sobol_sequence_indices[cell_ix] = variability_index;
}
