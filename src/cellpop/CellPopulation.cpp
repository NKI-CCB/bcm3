#include "Utils.h"
#include "Cell.h"
#include "CellPopulation.h"
#include "Experiment.h"
#include "SBMLModel.h"
#include "VariabilityPseudoRandomIterator.h"

CellPopulation::CellPopulation()
	: active_cells(0)
{
}

CellPopulation::~CellPopulation()
{
}

void CellPopulation::Allocate(const Experiment* experiment, SBMLModel* model, size_t max_cells, std::string solver_type)
{
	cells.resize(max_cells);
	cell_parents.resize(max_cells);
	for (size_t i = 0; i < max_cells; i++) {
		cells[i] = new Cell(model, experiment);
		cells[i]->AllocateSolver(solver_type);
		cell_parents[i] = std::numeric_limits<size_t>::max();
	}
}

void CellPopulation::Reset()
{
	active_cells = 0;
	for (size_t i = 0; i < cell_parents.size(); i++) {
		cell_parents[i] = std::numeric_limits<size_t>::max();
	}
}

size_t CellPopulation::AddNewCell(VariabilityPseudoRandomIterator* variability_iterator, double time, Cell* parent, int child_ix, const VectorReal& transformed_variables, bool entry_time_variable, bool request_synchronization, size_t initial_number_of_cells)
{
	if (active_cells == cells.size()) {
		return std::numeric_limits<size_t>::max();
	}

	bool result = true;

	size_t new_cell_ix = active_cells;
	Cell* cell = cells[new_cell_ix];
	active_cells++;

	// Store parent hierarchy
	if (parent) {
		size_t parent_ix = std::find(cells.begin(), cells.end(), parent) - cells.begin();
		cell_parents[new_cell_ix] = parent_ix;
	}

	// Set initial conditions for the new cell
	if (parent) {
		result &= cell->SetInitialConditionsFromOtherCell(parent);
	} else {
		// TODO
		//result &= cell->SetInitialConditionsFromModel(experiment->set_species_map, experiment->set_init_map, experiment->ratio_active_map, experiment->ratio_inactive_map, experiment->ratio_total_active, experiment->ratio_total_inactive, experiment->transformed_sampled_parameters, time);
		result &= cell->SetInitialConditionsFromModel();
	}

	size_t sobol_sequence_ix = 0;
	const VectorReal* sobol_values = nullptr;
	if (variability_iterator) {
		if (parent) {
			size_t parent_ix = cell_parents[new_cell_ix];

			int generation = 0;
			size_t grandparent_ix = cell_parents[parent_ix];
			while (1) {
				if (grandparent_ix == std::numeric_limits<size_t>::max()) {
					break;
				} else {
					grandparent_ix = cell_parents[grandparent_ix];
					generation++;
				}
			}

			sobol_sequence_ix = variability_iterator->GetIndex(parent_ix) * 2 + child_ix;
			while (generation > 0) {
				sobol_sequence_ix += initial_number_of_cells * (1 << generation);
				generation--;
			}
			if (sobol_sequence_ix >= variability_iterator->GetMaxIndex()) {
				return std::numeric_limits<size_t>::max();
			}
		} else {
			sobol_sequence_ix = new_cell_ix;
		}

		variability_iterator->SetIndex(new_cell_ix, sobol_sequence_ix);
		sobol_values = &variability_iterator->GetValue(sobol_sequence_ix);
	}

	// Initialize the cell for simulation
	result &= cell->Initialize(time, transformed_variables, sobol_values, entry_time_variable, request_synchronization);

	if (!result) {
		return std::numeric_limits<size_t>::max();
	} else {
		return new_cell_ix;
	}
}

size_t CellPopulation::CountCellsAtTime(double time, ESynchronizeCellTrajectory synchronize, bool count_only_mitotic) const
{
	size_t count = 0;
	for (size_t i = 0; i < active_cells; i++) {
		if (cells[i]->CellAliveAtTime(time, synchronize)) {
			if (count_only_mitotic) {
				if (cells[i]->EnteredMitosis()) {
					count++;
				}
			} else {
				count++;
			}
		}
	}
	return count;
}
