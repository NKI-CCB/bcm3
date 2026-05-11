#include "Utils.h"
#include "Cell.h"
#include "CellPopulation.h"
#include "DataLikelihoodBase.h"
#include "Experiment.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"
#include "SolverCodeGenerator.h"
#include "Timer.h"
#include "VariabilityDescription.h"
#include "VariabilityPseudoRandomIterator.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <fstream>

#define DUMP_CELL_TIMINGS 0

static const int output_trajectory_num_timepoints = 500;

Experiment::SimulationTimepoints::SimulationTimepoints()
	: data_likelihood_ix(std::numeric_limits<size_t>::max())
	, time(std::numeric_limits<Real>::quiet_NaN())
	, time_ix(std::numeric_limits<size_t>::max())
	, species_ix(std::numeric_limits<size_t>::max())
	, synchronize(ESynchronizeCellTrajectory::None)
{
}

void experiment_evaluation_worker(Experiment* experiment, size_t threadIndex)
{
	experiment->AuxWorkerFunction(threadIndex);
}

Experiment::AuxEvaluation::AuxEvaluation()
	: finished(false)
	, result(false)
	, terminate(false)
{
}

Experiment::Experiment(std::shared_ptr<const bcm3::VariableSet> varset, size_t evaluation_threads, bool store_simulation)
	: varset(varset)
	, cell_model(0)
	, initial_number_of_cells(0)
	, max_number_of_cells(0)
	, divide_cells(true)
	, entry_time_varix(std::numeric_limits<size_t>::max())
	, entry_time_nonsampled_varix(std::numeric_limits<size_t>::max())
	, fixed_entry_time(0.0)
	, synchronization_time_offset_varix(std::numeric_limits<size_t>::max())
	, fixed_synchronization_time_offset(0.0)
	, trailing_simulation_time(0.0)
	, simulate_past_chromatid_separation_time(0.0)
	, solver_min_timestep(4 * std::numeric_limits<float>::epsilon())
	, solver_max_timestep(std::numeric_limits<Real>::infinity())
	, solver_max_steps(10000)
	, solver_abs_tol(4 * std::numeric_limits<float>::epsilon())
	, solver_rel_tol(4 * std::numeric_limits<float>::epsilon())
	, derivative(nullptr)
	, jacobian(nullptr)
	, requested_duration(EPhaseDuration::None)
	, any_requested_synchronization(false)
	, store_simulation(store_simulation)
	, cell_submit_count(0)
	, cell_done_count(0)
	, aux_target_time(std::numeric_limits<Real>::quiet_NaN())
	, min_start_time(std::numeric_limits<Real>::quiet_NaN())
	, max_achieved_time(std::numeric_limits<Real>::quiet_NaN())
{
	if (evaluation_threads > 1) {
		AuxEvaluationThreads.resize(evaluation_threads);
		for (size_t i = 0; i < AuxEvaluationThreads.size(); i++) {
			AuxEvaluationThreads[i] = new AuxEvaluation;
			AuxEvaluationThreads[i]->thread = std::make_shared<std::thread>(experiment_evaluation_worker, this, i);
		}
	}
}

Experiment::~Experiment()
{
	for (size_t i = 0; i < AuxEvaluationThreads.size(); i++) {
		AuxEvaluation& e = *AuxEvaluationThreads[i];
		{
			std::lock_guard<std::mutex> lock(e.mutex);
			e.finished = false;
			e.terminate = true;
		}
		e.condition.notify_one();
		e.thread->join();
		delete AuxEvaluationThreads[i];
	}
}

bool Experiment::AddSimulationTimepoints(DataLikelihoodBase* dl, Real time, size_t time_ix, size_t species_ix, ESynchronizeCellTrajectory synchronize)
{
	SimulationTimepoints st;

	bool found = false;
	st.data_likelihood_ix = 0;
	for (auto dli = data_likelihoods.begin(); dli != data_likelihoods.end(); ++dli) {
		if (dli->get() == dl) {
			found = true;
			break;
		} else {
			st.data_likelihood_ix++;
		}
	}
	if (!found) {
		return false;
	}

	st.time = time;
	st.time_ix = time_ix;
	st.species_ix = species_ix;
	st.synchronize = synchronize;
	simulation_timepoints.push_back(st);

	if (synchronize != ESynchronizeCellTrajectory::None) {
		any_requested_synchronization = true;
	}

	return true;
}

bool Experiment::AddRequestedDuration(DataLikelihoodBase* dl, EPhaseDuration duration)
{
	requested_duration = duration;
	return true;
}

bool Experiment::AddNonSampledParameters(const std::vector<std::string>& variable_names)
{
	non_sampled_parameter_names = variable_names;
	non_sampled_parameters.setConstant(variable_names.size(), std::numeric_limits<Real>::quiet_NaN());
	return cell_model.AddNonSampledParameters(variable_names);
}

void Experiment::SetNonSampledParameters(const VectorReal& values)
{
	ASSERT(values.size() == non_sampled_parameter_names.size());
	non_sampled_parameters = values;
}

bool Experiment::PostInitialize()
{
	if (!cell_model.Initialize()) {
		return false;
	}
	for (size_t i = 0; i < data_likelihoods.size(); i++) {
		data_likelihoods[i]->PostInitialize(*varset.get(), non_sampled_parameter_names);
	}
	for (size_t i = 0; i < cell_variabilities.size(); i++) {
		cell_variabilities[i]->PostInitialize(varset, non_sampled_parameter_names, cell_model);
	}

	entry_time_varix = varset->GetVariableIndex(entry_time_varname, false);
	if (entry_time_varix == std::numeric_limits<size_t>::max()) {
		auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), entry_time_varname);
		if (it != non_sampled_parameter_names.end()) {
			entry_time_nonsampled_varix = it - non_sampled_parameter_names.begin();
		} else {
			try {
				fixed_entry_time = boost::lexical_cast<Real>(entry_time_varname);
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for entry_time \"%s\", and could also not cast it to a constant real value: %s", entry_time_varname.c_str(), e.what());
				return false;
			}
		}
	}

	if (!synchronization_time_offset_varname.empty()) {
		synchronization_time_offset_varix = varset->GetVariableIndex(synchronization_time_offset_varname, false);
		if (synchronization_time_offset_varix == std::numeric_limits<size_t>::max()) {
			try {
				fixed_entry_time = boost::lexical_cast<Real>(synchronization_time_offset_varname);
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Synchronization time offset was specified as \"%s\", but could not find variable and also could not cast it to a constant real value: %s", synchronization_time_offset_varname.c_str(), e.what());
				return false;
			}
		}
	} else {
		synchronization_time_offset_varix = std::numeric_limits<size_t>::max();
		fixed_synchronization_time_offset = 0.0;
	}

	Real simulation_begin_time = 0;
	Real simulation_end_time = 0;

	std::vector<Real> request_simulation_timepoints;
	for (size_t i = 0; i < simulation_timepoints.size(); i++) {
		request_simulation_timepoints.push_back(simulation_timepoints[i].time);
		simulation_end_time = std::max(simulation_end_time, simulation_timepoints[i].time);
	}

	if (store_simulation) {
		// Allocate space to store the trajectories
		output_trajectories_timepoints.resize(output_trajectory_num_timepoints);
		for (size_t i = 0; i < output_trajectories_timepoints.size(); i++) {
			output_trajectories_timepoints(i) = simulation_begin_time + (simulation_end_time - simulation_begin_time) * i / (Real)(output_trajectory_num_timepoints - 1);
			request_simulation_timepoints.push_back(output_trajectories_timepoints(i));
		}

		simulated_trajectories.resize(max_number_of_cells);
		for (size_t i = 0; i < simulated_trajectories.size(); i++) {
			simulated_trajectories[i].resize(output_trajectory_num_timepoints, VectorReal::Constant(cell_model.GetNumSimulatedSpecies(), std::numeric_limits<Real>::quiet_NaN()));
		}
	}

	std::sort(request_simulation_timepoints.begin(), request_simulation_timepoints.end());

	simulation_timepoints_vector = VectorReal(request_simulation_timepoints.size());
	for (size_t i = 0; i < request_simulation_timepoints.size(); i++) {
		simulation_timepoints_vector(i) = request_simulation_timepoints[i];
	}

	if (Cell::use_generated_code) {
		// Generate code for the model/variable combination
		code_generator = std::make_unique<SolverCodeGenerator>();

		std::string codegen_name = model_filename;
		codegen_name.at(codegen_name.find_last_of('.')) = '_';
		codegen_name += "_" + varset->GetVarsetName() + "_" + Name + "_codegen/";
		if (code_generator->GenerateAndCompileSolverCode(this, cell_model, codegen_name)) {
			derivative = code_generator->GetGeneratedDerivative();
			jacobian = code_generator->GetGeneratedJacobian();
		} else {
			return false;
		}
	}

	// Allocate space for cells
	population = std::make_unique<CellPopulation>();
	population->Allocate(this, &cell_model, max_number_of_cells, solver_type);

	return true;
}

bool Experiment::EvaluateLogProbability(size_t threadix, const VectorReal& values, const VectorReal& transformed_values, Real& logp)
{
	for (size_t i = 0; i < data_likelihoods.size(); i++) {
		data_likelihoods[i]->Reset();
	}
	if (store_simulation) {
		for (size_t i = 0; i < simulated_trajectories.size(); i++) {
			for (size_t j = 0; j < simulated_trajectories[i].size(); j++) {
				simulated_trajectories[i][j].setConstant(std::numeric_limits<Real>::quiet_NaN());
			}
		}
	}

	population->Reset();
	max_achieved_time = -std::numeric_limits<Real>::infinity();
	bool result = Simulate(transformed_values);

	if (result) {
		// Report durations
		if (requested_duration != EPhaseDuration::None) {
			for (size_t i = 0; i < population->GetActiveCount(); i++) {
				Real duration = population->GetCell(i)->GetDuration(requested_duration);
				for (size_t j = 0; j < data_likelihoods.size(); j++) {
					data_likelihoods[j]->NotifyDuration(i, duration);
				}
			}
		}

		// Synchronize cells if necessary
		Real time_offset = 0.0;
		if (synchronization_time_offset_varix != std::numeric_limits<size_t>::max()) {
			time_offset = transformed_sampled_parameters[synchronization_time_offset_varix];
		} else {
			time_offset = fixed_synchronization_time_offset;
		}

		// Report data references
		if (any_requested_synchronization) {
			for (int synchronize_i = 0; synchronize_i <= (int)ESynchronizeCellTrajectory::None; synchronize_i++) {
				// Need to reset the timepoint iteration for each different synchronization timepoint
				for (size_t i = 0; i < population->GetActiveCount(); i++) {
					population->GetCell(i)->RestartInterpolationIteration();
				}
				for (size_t ti = 0; ti < simulation_timepoints.size(); ti++) {
					const SimulationTimepoints& st = simulation_timepoints[ti];
					if (st.species_ix != std::numeric_limits<size_t>::max() && (int)st.synchronize == synchronize_i) {
						size_t population_size = population->CountCellsAtTime(st.time + time_offset, st.synchronize, false);
						size_t mitotic_population_size = population->CountCellsAtTime(st.time + time_offset, st.synchronize, true);
						for (size_t i = 0; i < population->GetActiveCount(); i++) {
							Cell* cell = population->GetCell(i);
							Real x = cell->GetInterpolatedSpeciesValue(st.time + time_offset, st.species_ix, st.synchronize);
							if (x == x) {
								data_likelihoods[st.data_likelihood_ix]->NotifySimulatedValue(st.time_ix, x, st.species_ix, i, population_size, mitotic_population_size, 0, cell->EnteredMitosis(), i >= initial_number_of_cells);
							}
						}
					}
				}
			}
		} else {
			for (size_t ti = 0; ti < simulation_timepoints.size(); ti++) {
				const SimulationTimepoints& st = simulation_timepoints[ti];
				if (st.species_ix != std::numeric_limits<size_t>::max()) {
					size_t population_size = population->CountCellsAtTime(st.time + time_offset, st.synchronize, false);
					size_t mitotic_population_size = population->CountCellsAtTime(st.time + time_offset, st.synchronize, true);
					for (size_t i = 0; i < population->GetActiveCount(); i++) {
						Cell* cell = population->GetCell(i);
						Real x = cell->GetInterpolatedSpeciesValue(st.time + time_offset, st.species_ix, st.synchronize);
						if (x == x) {
							data_likelihoods[st.data_likelihood_ix]->NotifySimulatedValue(st.time_ix, x, st.species_ix, i, population_size, mitotic_population_size, 0, cell->EnteredMitosis(), i >= initial_number_of_cells);
						}
					}
				}
			}
		}

		if (store_simulation) {
			for (size_t i = 0; i < population->GetActiveCount(); i++) {
				population->GetCell(i)->RestartInterpolationIteration();
			}

			Real simulation_begin_time = min_start_time;
			Real simulation_end_time = max_achieved_time;

			for (size_t i = 0; i < output_trajectories_timepoints.size(); i++) {
				output_trajectories_timepoints(i) = simulation_begin_time + (simulation_end_time - simulation_begin_time) * i / (Real)(output_trajectory_num_timepoints - 1);
			}

			for (size_t ti = 0; ti < output_trajectories_timepoints.size(); ti++) {
				Real output_t = output_trajectories_timepoints(ti);
				for (size_t i = 0; i < population->GetActiveCount(); i++) {
					Cell* cell = population->GetCell(i);
					for (size_t j = 0; j < cell_model.GetNumODEIntegratedSpecies(); j++) {
						size_t simix = cell_model.GetSimulatedSpeciesFromODEIntegratedSpecies(j);
						simulated_trajectories[i][ti](simix) = cell->GetInterpolatedSpeciesValue(output_t, j, ESynchronizeCellTrajectory::None);
					}
					for (size_t j = 0; j < cell_model.GetNumConstantSpecies(); j++) {
						size_t simix = cell_model.GetSimulatedSpeciesFromConstantSpecies(j);
						simulated_trajectories[i][ti](simix) = cell->GetInterpolatedSpeciesValue(output_t, cell_model.GetNumODEIntegratedSpecies() + j, ESynchronizeCellTrajectory::None);
					}
					for (size_t j = 0; j < treatment_trajectories.size(); j++) {
						size_t simix = cell_model.GetSimulatedSpeciesByName(cell_model.GetConstantSpeciesName(treatment_trajectories_species_ix[j]));
						simulated_trajectories[i][ti](simix) = treatment_trajectories[j]->GetConcentration(output_t, cell->GetCreationTime());
					}
				}
			}
		}

		// Evaluate the requested time points
		logp = 0.0;
		for (size_t i = 0; i < data_likelihoods.size(); i++) {
			Real dl_logp = 0.0;
			if (!data_likelihoods[i]->Evaluate(values, transformed_values, non_sampled_parameters, dl_logp)) {
				result = false;
				break;
			}
			logp += dl_logp;
		}
	} else {
		logp = -std::numeric_limits<Real>::infinity();
	}

	Real max_solver_steps = 0;
	Real min_solver_step_size = std::numeric_limits<Real>::infinity();

	for (size_t i = 0; i < population->GetActiveCount(); i++) {
		Cell* cell = population->GetCell(i);
		max_solver_steps = (std::max)(max_solver_steps, (Real)cell->GetSolverSteps());
		min_solver_step_size = (std::min)(min_solver_step_size, (Real)cell->GetSolverMinStepSize());
	}

	solver_stats.push_back(std::tuple<Real, Real, bool>(max_solver_steps, min_solver_step_size, result));

	return true;
}

void Experiment::DumpCVodeStatistics(const std::string& output_folder)
{
	std::string fn = output_folder + std::string("/cvode_statistics_") + Name + std::string(".tsv");
	FILE* f = fopen(fn.c_str(), "w");
	fprintf(f, "cvode_steps\tmin_step_size\tsucceeded\n");
	for (size_t i = 0; i < solver_stats.size(); i++) {
		fprintf(f, "%g\t%g\t%d\n", std::get<0>(solver_stats[i]), std::get<1>(solver_stats[i]), std::get<2>(solver_stats[i]));
	}
	fclose(f);
}

std::unique_ptr<Experiment> Experiment::Create(const boost::property_tree::ptree& xml_node, std::shared_ptr<const bcm3::VariableSet> varset, const boost::program_options::variables_map& vm, bcm3::RNG& rng, size_t evaluation_threads, bool store_simulation)
{
	std::unique_ptr<Experiment> experiment;

	try {
		experiment = std::make_unique<Experiment>(varset, evaluation_threads, store_simulation);
		if (!experiment->Load(xml_node, vm)) {
			experiment.reset();
		} else {
			experiment->rng.Seed(rng.GetUnsignedInt());
		}
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		experiment.reset();
	}

	return experiment;
}

bool Experiment::Load(const boost::property_tree::ptree& xml_node, const boost::program_options::variables_map& vm)
{
	Name = xml_node.get<std::string>("<xmlattr>.name");

	model_filename = xml_node.get<std::string>("<xmlattr>.model_file");
	std::string data_file = xml_node.get<std::string>("<xmlattr>.data_file", "");

	solver_type = xml_node.get<std::string>("<xmlattr>.solver_type", "CVODE");
	solver_min_timestep = xml_node.get<Real>("<xmlattr>.solver_min_timestep", 1e-8);
	solver_max_timestep = xml_node.get<Real>("<xmlattr>.solver_max_timestep", std::numeric_limits<Real>::infinity());
	solver_max_steps = xml_node.get<size_t>("<xmlattr>.solver_max_steps", 10000);
	solver_abs_tol = xml_node.get<Real>("<xmlattr>.solver_absolute_tolerance", 4 * std::numeric_limits<float>::epsilon());
	solver_rel_tol = xml_node.get<Real>("<xmlattr>.solver_relative_tolerance", 4 * std::numeric_limits<float>::epsilon());

	// Load & initialize the model
	if (!cell_model.LoadSBML(model_filename)) {
		return false;
	}
	if (!cell_model.SetVariableSet(varset)) {
		return false;
	}

	//if any ratios are present they are mapped here
	// if any initial values are present as parameters they are mapped here to the ODE model
	size_t num_variables = varset->GetNumVariables();
	std::vector<std::string> variable_names = varset->GetAllVariableNames();

	for(int i = 0; i < num_variables; i++){
		std::string name_var = variable_names[i];
		if(name_var.substr(0,8).compare("species_") == 0){
			set_init_map[GetODEIntegratedSpeciesByName(name_var.substr(8), true)] = i;
		}

		if(name_var.substr(0,6).compare("ratio_") == 0){
			// Check that the total variable is also present
			bool also_total_var = false;

			for(int j = 0; j < num_variables; j++){
				std::string name_var_internal = variable_names[j];
				if(name_var_internal.substr(0,6).compare("total_") == 0 && name_var.substr(6).compare(name_var_internal.substr(6)) == 0){
					also_total_var = true;
					std::vector<int> v = {i, j};
					size_t active_species = GetODEIntegratedSpeciesByName("active_" + name_var.substr(6), true);
					size_t inactive_species = GetODEIntegratedSpeciesByName("inactive_" + name_var.substr(6), true);

					if(active_species == std::numeric_limits<size_t>::max()){
						LOG("The ratio variable \"%s\" has been specified; \"ratio_\" variables are used to describe ratios of initial conditions for two species; but the corresponding active species \"active_%s\" is not found in the model.", name_var.c_str(), name_var.substr(6).c_str());
						return false;
					}

					if(inactive_species == std::numeric_limits<size_t>::max()){
						LOG("The ratio variable \"%s\" has been specified; \"ratio_\" variables are used to describe ratios of initial conditions for two species; but the corresponding inactive species \"inactive_%s\" is not found in the model.", name_var.c_str(), name_var.substr(6).c_str());
						return false;
					}

					ratio_active_map[active_species] = v;
					ratio_inactive_map[inactive_species] = v;
				}
			}

			if (!also_total_var) {
				size_t active_species = GetODEIntegratedSpeciesByName("active_" + name_var.substr(6), true);
				size_t inactive_species = GetODEIntegratedSpeciesByName("inactive_" + name_var.substr(6), true);

				if(active_species == std::numeric_limits<size_t>::max()){
					LOG("The ratio variable \"%s\" has been specified; \"ratio_\" variables are used to describe ratios of initial conditions for two species; but the corresponding active species \"active_%s\" is not found in the model.", name_var.c_str(), name_var.substr(6).c_str());
					return false;
				}

				if(inactive_species == std::numeric_limits<size_t>::max()){
					LOG("The ratio variable \"%s\" has been specified; \"ratio_\" variables are used to describe ratios of initial conditions for two species; but the corresponding inactive species \"inactive_%s\" is not found in the model.", name_var.c_str(), name_var.substr(6).c_str());
					return false;
				}

				std::vector<size_t> v = {(size_t)i, active_species, inactive_species};

				ratio_total_active[active_species] = v;
				ratio_total_inactive[inactive_species] = v;
			}
		}
	}

	initial_number_of_cells = xml_node.get<size_t>("<xmlattr>.num_cells", 1);
	max_number_of_cells = xml_node.get<size_t>("<xmlattr>.max_cells", 20);
	divide_cells = xml_node.get<bool>("<xmlattr>.divide_cells", true);
	trailing_simulation_time = xml_node.get<Real>("<xmlattr>.trailing_simulation_time", 0.0);
	simulate_past_chromatid_separation_time = xml_node.get<Real>("<xmlattr>.simulate_past_chromatid_separation_time", 0.0);

	// Load the species/parameter setting
	try {
		BOOST_FOREACH(const boost::property_tree::ptree::value_type& var, xml_node.get_child("")) {
			if (var.first == "set_species") {
				std::string species_name = var.second.get<std::string>("<xmlattr>.species_name");
				size_t species_ix = cell_model.GetSimulatedSpeciesByName(species_name);
				if (species_ix == std::numeric_limits<size_t>::max()) {
					return false;
				}
				SetSpecies ss;
				ss.value = var.second.get<Real>("<xmlattr>.value");
				ss.begin_time = var.second.get<Real>("<xmlattr>.begin_time", -std::numeric_limits<Real>::infinity());
				ss.end_time = var.second.get<Real>("<xmlattr>.end_time", std::numeric_limits<Real>::infinity());
				ss.old_value = 0;//TODO; read from SBML file
				set_species_map[species_ix] = ss;
			}
			if (var.first == "set_parameter") {
				std::string parameter_name = var.second.get<std::string>("<xmlattr>.parameter_name");
				Real value = var.second.get<Real>("<xmlattr>.value");
				if (!cell_model.AddFixedParameter(parameter_name, value)) {
					return false;
				}
			}
			if (var.first == "experiment_specific_parameter") {
				std::string parameter_name = var.second.get<std::string>("<xmlattr>.parameter_name");
				size_t parameter_ix = varset->GetVariableIndex(parameter_name);
				if (parameter_ix == std::numeric_limits<size_t>::max()) {
					return false;
				}
				std::string replacement_parameter_name = var.second.get<std::string>("<xmlattr>.replacement_parameter_name");
				size_t replacement_parameter_ix = varset->GetVariableIndex(replacement_parameter_name);
				if (replacement_parameter_ix == std::numeric_limits<size_t>::max()) {
					return false;
				}
				experiment_specific_parameter_map[parameter_ix] = replacement_parameter_ix;
			}
		}
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	// Load the data variability specifications
	any_requested_synchronization = false;
	try {
		BOOST_FOREACH(const boost::property_tree::ptree::value_type & var, xml_node.get_child("")) {
			if (var.first == "cell_variability") {
				std::unique_ptr<VariabilityDescription> vd = VariabilityDescription::Create(var.second);
				if (vd) {
					cell_variabilities.push_back(std::move(vd));
				} else {
					return false;
				}
			}
		}
	} catch (boost::property_tree::ptree_error& e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	// Load the data file
	if (!data_file.empty()) {
		bcm3::NetCDFDataFile file;
		if (!file.Open(data_file, false)) {
			return false;
		}

		// Load the data references and variability specifications
		any_requested_synchronization = false;
		try {
			BOOST_FOREACH(const boost::property_tree::ptree::value_type& var, xml_node.get_child("")) {
				if (var.first == "data") {
					std::unique_ptr<DataLikelihoodBase> dl = DataLikelihoodBase::Create(var.second, AuxEvaluationThreads.size());
					data_likelihoods.push_back(std::move(dl)); // Need to have the data likelihood in the member list already
					if (!(*data_likelihoods.rbegin())->Load(var.second, this, *varset.get(), file, vm)) {
						data_likelihoods.pop_back();
						return false;
					}
				}
				if (var.first == "treatment_trajectory") {
					std::string species_name = var.second.get<std::string>("<xmlattr>.species_name");
					size_t constant_species_ix = cell_model.GetConstantSpeciesByName(species_name);
					if (constant_species_ix == std::numeric_limits<size_t>::max()) {
						LOGERROR("Cannot find \"%s\" as a constant species for treatment trajectory (the species needs to be constant).", species_name.c_str());
						return false;
					}
					treatment_trajectories_species_ix.push_back(constant_species_ix);

					std::unique_ptr<TreatmentTrajectory> traj = TreatmentTrajectory::Create(var.second);
					if (traj->Load(var.second, this, file)) {
						treatment_trajectories.push_back(std::move(traj));
					} else {
						return false;
					}
				}
			}
		} catch (boost::property_tree::ptree_error &e) {
			LOGERROR("Error parsing likelihood file: %s", e.what());
			return false;
		}

		// Sort the requested timepoints by time
		for (size_t i = 1; i < simulation_timepoints.size(); i++) {
			for (size_t j = 0; j < simulation_timepoints.size() - i; j++) {
				if (simulation_timepoints[j].time > simulation_timepoints[j + 1].time) {
					std::swap(simulation_timepoints[j], simulation_timepoints[j + 1]);
				}
			}
		}
	} else {
		LOG("No data file specified - adding fixed timepoints to simulate model");
		for (int i = 0; i < 2; i++) {
			SimulationTimepoints st;
			st.data_likelihood_ix = std::numeric_limits<size_t>::max();
			st.species_ix = std::numeric_limits<size_t>::max();
			st.time = i * 2000.0;
			st.time_ix = i;
			simulation_timepoints.push_back(st);
		}
	}

	return Initialize(xml_node);
}

bool Experiment::Initialize(const boost::property_tree::ptree& xml_node)
{
	entry_time_varname = xml_node.get<std::string>("<xmlattr>.entry_time");
	synchronization_time_offset_varname = xml_node.get<std::string>("<xmlattr>.synchronization_time_offset", "");

	transformed_sampled_parameters.setConstant(varset->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());

	if (cell_variabilities.size() > 0) {
		size_t total_dimensions = 0;
		for (size_t i = 0; i < cell_variabilities.size(); i++) {
			total_dimensions += cell_variabilities[i]->GetNumDimensions();
		}
		variability_iterator = std::make_unique<VariabilityPseudoRandomIterator>();
		variability_iterator->Initialize(total_dimensions, initial_number_of_cells, max_number_of_cells);
	}

	return true;
}

bool Experiment::Simulate(const VectorReal& transformed_values)
{
	for (size_t i = 0; i < varset->GetNumVariables(); i++) {
		transformed_sampled_parameters[i] = transformed_values[i];
	}
	for (std::map<size_t, size_t>::iterator espi = experiment_specific_parameter_map.begin(); espi != experiment_specific_parameter_map.end(); ++espi) {
		transformed_sampled_parameters[espi->first] = transformed_sampled_parameters[espi->second];
	}

	Real entry_time;
	if (entry_time_varix != std::numeric_limits<size_t>::max()) {
		entry_time = transformed_sampled_parameters[entry_time_varix];
	} else if (entry_time_nonsampled_varix != std::numeric_limits<size_t>::max()) {
		entry_time = non_sampled_parameters[entry_time_nonsampled_varix];
	} else {
		entry_time = fixed_entry_time;
	}

	const Real dt = 120.0;
	Real simulation_end_time;
	if (!simulation_timepoints.empty()) {
		simulation_end_time = simulation_timepoints.rbegin()->time + trailing_simulation_time;
	} else {
		simulation_end_time = trailing_simulation_time;
	}

	min_start_time = std::numeric_limits<Real>::infinity();
	if (initial_number_of_cells > 1) {
		for (size_t i = 0; i < initial_number_of_cells; i++) {
			Real cell_start_time = entry_time;
			population->AddNewCell(variability_iterator.get(), cell_start_time, NULL, -1, transformed_sampled_parameters, true, any_requested_synchronization, initial_number_of_cells);
			min_start_time = std::min(cell_start_time, min_start_time);
		}
	} else {
		population->AddNewCell(variability_iterator.get(), entry_time, NULL, -1, transformed_sampled_parameters, false, any_requested_synchronization, initial_number_of_cells);
		min_start_time = entry_time;
	}

	if (min_start_time < -7.0 * 24.0 * 60.0 * 60.0) {
		//printf("Minimum start time too long: %g", minimum_start_time);
		return false;
	}

	for (size_t i = 0; i < initial_number_of_cells; i++) {
		for (auto dli = data_likelihoods.begin(); dli != data_likelihoods.end(); ++dli) {
			(*dli)->NotifyStartingCells(i);
		}
	}

	if (!ParallelSimulation(simulation_end_time)) {
		return false;
	}

	return true;
}

bool Experiment::ParallelSimulation(Real target_time)
{
	bool result = true;
	Real achieved_time = 0.0;

	if (!AuxEvaluationThreads.empty()) {
		cells_to_process_lock.lock();
		aux_target_time = target_time;
		for (size_t i = 0; i < population->GetActiveCount(); i++) {
			cells_to_process.push(i);
		}
		cells_to_process_lock.unlock();

		{
			std::unique_lock<std::mutex> lock(all_done_mutex);
			cell_submit_count = population->GetActiveCount();
			cell_done_count = 0;
		}

		StartAuxThreads();
		result = WaitAuxThreads();
	} else {
		for (size_t i = 0; i < population->GetActiveCount(); i++) {
			result &= SimulateCell(i, target_time, achieved_time, 0);
			if (!result) {
				break;
			} else {
				max_achieved_time = std::max(achieved_time, max_achieved_time);
			}
		}
	}

	return result;
}

bool Experiment::SimulateCell(size_t i, Real target_time, Real& achieved_time, size_t eval_thread)
{
	bool die = false, divide = false;
	achieved_time = 0.0;
	if (!population->GetCell(i)->Simulate(target_time, simulate_past_chromatid_separation_time, simulation_timepoints_vector, die, divide, achieved_time)) {
		return false;
	}

	if (divide_cells && divide && achieved_time < target_time) {
		size_t cell1, cell2;
		// This is the only place where the number of active cells is increased; the readers don't need to lock it
		{
			std::lock_guard<std::mutex> lock(cell_vector_mutex);

			// Create two new cells
			Cell* parent = population->GetCell(i);
			cell1 = population->AddNewCell(variability_iterator.get(), achieved_time, parent, 0, transformed_sampled_parameters, false, any_requested_synchronization, initial_number_of_cells);
			cell2 = population->AddNewCell(variability_iterator.get(), achieved_time, parent, 1, transformed_sampled_parameters, false, any_requested_synchronization, initial_number_of_cells);
			if (cell1 == std::numeric_limits<size_t>::max() || cell2 == std::numeric_limits<size_t>::max()) {
				return false;
			}
			for (auto dli = data_likelihoods.begin(); dli != data_likelihoods.end(); ++dli) {
				(*dli)->NotifyParents(i, cell1);
				(*dli)->NotifyParents(i, cell2);
			}
		}

		{
			std::unique_lock<std::mutex> lock(all_done_mutex);
			cell_submit_count += 2;

			cells_to_process_lock.lock();
			cells_to_process.push(cell1);
			cells_to_process.push(cell2);
			cells_to_process_lock.unlock();
		}

		// Try to wake up one other aux threads if any is sleeping
		for (size_t i = 0; i < AuxEvaluationThreads.size(); i++) {
			if (i != eval_thread) {
				AuxEvaluation& e = *AuxEvaluationThreads[i];
				if (e.mutex.try_lock()) {
					if (e.finished) {
						e.finished = false;
						e.mutex.unlock();
						e.condition.notify_one();
						break;
					} else {
						e.mutex.unlock();
					}
				}
			}
		}
	}

	return true;
}

void Experiment::StartAuxThreads()
{
	for (size_t i = 0; i < AuxEvaluationThreads.size(); i++) {
		AuxEvaluation& e = *AuxEvaluationThreads[i];
		{
			std::lock_guard<std::mutex> lock(e.mutex);
			e.result = true;
			e.finished = false;
		}
		e.condition.notify_one();
	}
}

bool Experiment::WaitAuxThreads()
{
	bool result = true;
	std::unique_lock<std::mutex> lock(all_done_mutex);
	while (cell_done_count != cell_submit_count) {
		all_done_condition.wait(lock);
	}
	for (size_t i = 0; i < AuxEvaluationThreads.size(); i++) {
		AuxEvaluation& e = *AuxEvaluationThreads[i];
		result &= e.result;
	}
	return result;
}

void Experiment::AuxWorkerFunction(size_t threadIndex)
{
	AuxEvaluation& e = *AuxEvaluationThreads[threadIndex];

	while (1) {
		{
			std::unique_lock<std::mutex> lock(e.mutex);
			while (e.finished) {
				e.condition.wait(lock);
			}
		}

		if (e.terminate) {
			break;
		}

		bool queue_empty = false;
		size_t cells_processed = 0;
		while (!queue_empty) {
			size_t cell_ix;
			cells_to_process_lock.lock();
			if (cells_to_process.empty()) {
				cells_to_process_lock.unlock();
				queue_empty = true;
			} else {
				cell_ix = cells_to_process.front();
				cells_to_process.pop();
				cells_to_process_lock.unlock();

#if DUMP_CELL_TIMINGS
				bcm3::Timer timer;
				timer.Start();
#endif

				Real achieved_time;
				e.result &= SimulateCell(cell_ix, aux_target_time, achieved_time, threadIndex);
				cells_processed++;

				cells_to_process_lock.lock();
				max_achieved_time = std::max(achieved_time, max_achieved_time);
				cells_to_process_lock.unlock();

#if DUMP_CELL_TIMINGS
				Real time = timer.GetElapsedSeconds();
				timings_lock.lock();
				(*cell_update_timings.rbegin())(cell_ix) = time;
				timings_lock.unlock();
#endif
			}


		}

		{
			std::lock_guard<std::mutex> lock(e.mutex);
			e.finished = true;
		}
		{
			std::unique_lock<std::mutex> lock(all_done_mutex);
			cell_done_count += cells_processed;
		}
		all_done_condition.notify_one();
	}
}
