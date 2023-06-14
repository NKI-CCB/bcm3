#pragma once

#include "CVODESolverDelay.h"
#include "Likelihood.h"
#include "RNG.h"
#include <boost/multi_array.hpp>

class SBMLModel;

class LikelihoodIncucytePopulation : public bcm3::Likelihood
{
public:
	LikelihoodIncucytePopulation(size_t sampling_threads, size_t evaluation_threads);
	~LikelihoodIncucytePopulation();

	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node, const boost::program_options::variables_map& vm);
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	MatrixReal GetSimulatedCellCount() const;
	MatrixReal GetSimulatedApoptoticCellCount() const;
	MatrixReal GetSimulatedDebris() const;
	MatrixReal GetSimulatedConfluence() const;
	MatrixReal GetSimulatedApoptosisMarker() const;
	VectorReal GetSimulatedCTB() const;

private:
	struct Well
	{
		Well() {}
		Well(size_t num_timepoints, size_t num_replicates);
		VectorReal cell_count;
		VectorReal apoptotic_cell_count;
		VectorReal debris;
		MatrixReal confluence;
		MatrixReal apoptosis_marker;
	};

	struct CellTreatment
	{
		std::vector<Well> drug_treatment;
		Well positive_control;
		Well negative_control;
		VectorReal ctb;
	};

	struct StochasticCell
	{
		Real cell_cycle_progression;
		Real size;
		Real apoptosis_progression;
		Real docetaxel_mitosos_block;
		bool debris;
	};

	struct ParallelData {
		Real proliferation_rate;
		Real apoptosis_rate;
		Real apoptosis_duration;
		Real apoptosis_remove_rate;

		Real cell_size;
		Real apoptotic_cell_size;
		Real debris_size;
		Real contact_inhibition_start;
		Real contact_inhibition_max_confluence;
		Real contact_inhibition_apoptosis_rate;

		Real drug_start_time;
		Real drug_effect_time;
		Real drug_proliferation_rate;
		Real drug_apoptosis_rate;
		VectorReal delay_y;

		bool pao;
		bool drug;

		bcm3::RNG rng;
	};

	struct Experiment {
		size_t experiment_ix;
		size_t num_timepoints;
		size_t num_replicates;
		size_t num_concentrations;

		VectorReal observed_timepoints;
		VectorReal concentrations;
		Real time_of_drug_treatment;
		Real time_of_ctb;
		Real seeding_density;

		CellTreatment observed_cell_treatment;
		std::vector< CellTreatment > modeled_cell_treatments;
	};

	bool SimulateWell(size_t threadix, Experiment& e, Well& w, bool pao, bool drug, bool single_drug, Real drug_proliferation_rate, Real drug_apoptosis_rate, const VectorReal& values);
	bool SimulateWellDeterministic(size_t threadix, Experiment& e, Well& w, bool pao, bool drug, bool single_drug, Real drug_proliferation_rate, Real drug_apoptosis_rate, const VectorReal& values);
	bool SimulateWellStochastic(size_t threadix, Experiment& e, Well& w, bool pao, bool drug, bool single_drug, Real drug_proliferation_rate, Real drug_apoptosis_rate, const VectorReal& values);
	bool CalculateDerivative(OdeReal t, const OdeReal* y, const std::vector< OdeReal >& history_t, const std::vector< OdeVectorReal >& history_y, size_t current_dci, OdeReal* dydt, void* user);

	std::shared_ptr<const bcm3::VariableSet> varset;
	size_t evaluation_threads;

	std::string drug_name;
	std::string cell_line;
	bool stochastic;
	bool use_pao_control;
	size_t num_combo_drugs;

	std::vector< Experiment > experiments;
	std::vector< CVODESolverDelay > solvers;
	std::vector< ParallelData > parallel_data;
};
