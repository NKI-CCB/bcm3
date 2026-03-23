#pragma once

class Cell;
class Experiment;
class SBMLModel;
class VariabilityPseudoRandomIterator;
enum class ESynchronizeCellTrajectory;

class CellPopulation
{
public:
	CellPopulation();
	~CellPopulation();

	void Allocate(const Experiment* experiment, SBMLModel* model, size_t max_cells, Real abs_tol, Real rel_tol);

	void Reset();
	size_t AddNewCell(VariabilityPseudoRandomIterator* variability_iterator, double time, Cell* parent, int child_ix, const VectorReal& transformed_variables, bool entry_time_variable, bool request_synchronization, size_t initial_number_of_cells);
	size_t CountCellsAtTime(double time, ESynchronizeCellTrajectory synchronize, bool count_only_mitotic) const;

	Cell* GetCell(size_t i) const { return cells[i]; }
	size_t GetCellParent(size_t i) const { return cell_parents[i]; }
	size_t GetActiveCount() const { return active_cells; }

private:
	std::vector<Cell*> cells;
	size_t active_cells;

	std::vector<size_t> cell_parents;
};
