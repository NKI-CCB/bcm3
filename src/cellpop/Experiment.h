#pragma once

#include <sundials/sundials_matrix.h>
#include <boost/program_options.hpp>
#include <condition_variable>

#include "CellPopulation.h"
#include "NetCDFDataFile.h"
#include "RNG.h"
#include "SBMLModel.h"
#include "Spinlock.h"
#include "TreatmentTrajectory.h"
#include "VariableSet.h"

extern "C" {
	typedef void (*derivative_fn)(OdeReal* out, const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
	typedef void (*jacobian_fn)(OdeMatrixReal& out, const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
}

class Cell;
class DataLikelihoodBase;
class VariabilityDescription;
class VariabilityPseudoRandomIterator;

enum class ESynchronizeCellTrajectory {
	DNAReplicationStart,
	PCNA_gfp_increase,
	NuclearEnvelopeBreakdown,
	AnaphaseOnset,
	None // None must be last
};

enum class EPhaseDuration {
	G1phase,
	Sphase,
	G2phase,
	NEBD_to_AnaphaseOnset,
	None // None must be last
};

class Experiment
{
public:
	Experiment(std::shared_ptr<const bcm3::VariableSet> varset, size_t evaluation_threads, bool store_simulation);
	virtual ~Experiment();

	bool AddSimulationTimepoints(DataLikelihoodBase* dl, Real time, size_t time_ix, size_t species_ix, ESynchronizeCellTrajectory synchronize);
	bool AddRequestedDuration(DataLikelihoodBase* dl, EPhaseDuration duration);
	bool AddNonSampledParameters(const std::vector<std::string>& variable_names);
	void SetNonSampledParameters(const VectorReal& values);
	bool PostInitialize();
	bool EvaluateLogProbability(size_t threadix, const VectorReal& values, const VectorReal& transformed_values, Real& logp);
	void DumpCVodeStatistics(const std::string& output_folder);

	static std::unique_ptr<Experiment> Create(const boost::property_tree::ptree& xml_node, std::shared_ptr<const bcm3::VariableSet> varset, const boost::program_options::variables_map& vm, bcm3::RNG& rng, size_t evaluation_threads, bool store_simulation);
	inline const std::string& GetName() const { return Name; }
	inline size_t GetNumODEIntegratedSpecies() const { return cell_model.GetNumODEIntegratedSpecies(); }
	inline size_t GetSimulatedSpeciesByName(const std::string& species_name) const { return cell_model.GetSimulatedSpeciesByName(species_name); }
	inline size_t GetODEIntegratedSpeciesByName(const std::string& species_name, bool log_error) const { return cell_model.GetODEIntegratedSpeciesByName(species_name, log_error); }
	inline size_t GetConstantSpeciesByName(const std::string& species_name, bool log_error) const { return cell_model.GetConstantSpeciesByName(species_name, log_error); }
	inline size_t GetMaxNumberOfCells() const { return max_number_of_cells; }
	inline const bcm3::VariableSet* GetVarset() const { return varset.get(); }

	inline size_t GetNumSpecies() const { return cell_model.GetNumSimulatedSpecies(); }
	inline std::string GetSpeciesName(size_t i) const { return cell_model.GetSimulatedSpeciesName(i); }
	inline size_t GetOutputNumCells() const { return std::min(population->GetActiveCount(), simulated_trajectories.size()); }
	inline const VectorReal& GetOutputTimepoints() const { return output_trajectories_timepoints; }
	inline const VectorReal& GetSimulatedTrajectory(size_t cell_i, size_t time_i) const { return simulated_trajectories[cell_i][time_i]; }
	inline size_t GetSimulatedCellParent(size_t cell_i) const { return population->GetCellParent(cell_i); }
	inline size_t GetNumData() const { return data_likelihoods.size(); }
	inline const DataLikelihoodBase* GetData(size_t i) const { return data_likelihoods[i].get(); }

protected:
	bool Load(const boost::property_tree::ptree& xml_node, const boost::program_options::variables_map& vm);
	bool GenerateAndCompileSolverCode(const std::string& codegen_name);
	bool Initialize(const boost::property_tree::ptree& xml_node);
	bool Simulate(const VectorReal& transformed_values);
	bool ParallelSimulation(Real target_time);
	bool SimulateCell(size_t i, Real target_time, Real& achieved_time, size_t eval_thread);

	void StartAuxThreads();
	bool WaitAuxThreads();
	void AuxWorkerFunction(size_t threadIndex);

	enum EInitialCellCycleDistribution {
		CCD_Uniform,
		CCD_Normal,
		CCD_Invalid,
	};

	std::string Name;

	// Settings
	std::string model_filename;
	SBMLModel cell_model;
	std::shared_ptr<const bcm3::VariableSet> varset;
	std::vector< std::unique_ptr<DataLikelihoodBase> > data_likelihoods;
	size_t initial_number_of_cells;
	size_t max_number_of_cells;
	bool divide_cells;
	std::string entry_time_varname;
	size_t entry_time_varix;
	size_t entry_time_nonsampled_varix;
	Real fixed_entry_time;
	std::string synchronization_time_offset_varname;
	size_t synchronization_time_offset_varix;
	Real fixed_synchronization_time_offset;
	Real trailing_simulation_time;
	Real simulate_past_chromatid_separation_time;
	std::vector< std::unique_ptr<VariabilityDescription> > cell_variabilities;

	std::vector< std::unique_ptr<TreatmentTrajectory> > treatment_trajectories;
	std::vector< size_t > treatment_trajectories_species_ix;

	struct SetSpecies {
		Real begin_time;
		Real end_time;
		Real value;
		Real old_value;
	};
	std::map<size_t, SetSpecies> set_species_map;
	std::map<size_t, size_t> set_init_map;
	std::map<size_t, std::vector<int>> ratio_active_map;
	std::map<size_t, std::vector<int>> ratio_inactive_map;
	std::map<size_t, std::vector<size_t>> ratio_total_active;
	std::map<size_t, std::vector<size_t>> ratio_total_inactive;
	std::map<size_t, size_t> experiment_specific_parameter_map;

	// Runtime variables
	std::vector<std::string> non_sampled_parameter_names;
	VectorReal transformed_sampled_parameters;
	VectorReal non_sampled_parameters;
	std::unique_ptr<CellPopulation> population;
	std::unique_ptr<VariabilityPseudoRandomIterator> variability_iterator;
	struct SimulationTimepoints {
		SimulationTimepoints();
		size_t data_likelihood_ix;
		Real time;
		size_t time_ix;
		size_t species_ix;
		ESynchronizeCellTrajectory synchronize;
	};
	std::vector<SimulationTimepoints> simulation_timepoints;
	VectorReal simulation_timepoints_vector;
	EPhaseDuration requested_duration;
	bool any_requested_synchronization;
	bool store_simulation;

	VectorReal output_trajectories_timepoints;
	std::vector<std::vector<VectorReal>> simulated_trajectories;

	Real abs_tol;
	Real rel_tol;

#if PLATFORM_WINDOWS
	HMODULE derivative_dll;
#else
	void* derivative_dll;
#endif
	derivative_fn derivative;
	jacobian_fn jacobian;
	mutable bcm3::RNG rng;

	std::vector< std::tuple<Real, Real, bool> > solver_stats;

	struct AuxEvaluation {
		AuxEvaluation();
		std::shared_ptr<std::thread> thread;
		std::mutex mutex;
		std::condition_variable condition;
		bool finished;
		bool result;
		bool terminate;
	};
	std::vector<AuxEvaluation*> AuxEvaluationThreads;
	Real aux_target_time;
	Real min_start_time;
	Real max_achieved_time;
	std::queue<size_t> cells_to_process;
	bcm3::spinlock cells_to_process_lock;
	std::vector<size_t> add_new_cell_parents;
	std::mutex cell_vector_mutex;
	size_t cell_submit_count;
	size_t cell_done_count;
	std::mutex all_done_mutex;
	std::condition_variable all_done_condition;

	bcm3::spinlock timings_lock;
	std::vector<VectorReal> cell_update_timings;

	friend class Cell;
	friend class Population;
	friend void experiment_evaluation_worker(Experiment* experiment, size_t threadIndex);
};
