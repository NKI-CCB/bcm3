#pragma once

#include <boost/random/sobol.hpp>

class VariabilityPseudoRandomIterator
{
public:
	VariabilityPseudoRandomIterator();
	~VariabilityPseudoRandomIterator();

	void Initialize(size_t dimensions, size_t initial_number_of_cells, size_t max_number_of_cells);

	const VectorReal& GetValue(size_t index) const { return sobol_sequence_values[index]; }

	size_t GetMaxIndex() const { return sobol_sequence_values.size(); }
	size_t GetIndex(size_t cell_ix) const { return sobol_sequence_indices[cell_ix]; }
	void SetIndex(size_t cell_ix, size_t variability_index);

private:
	std::shared_ptr<boost::random::sobol> sobol_sequence;
	std::vector<size_t> sobol_sequence_indices;
	std::vector<VectorReal> sobol_sequence_values;
};
