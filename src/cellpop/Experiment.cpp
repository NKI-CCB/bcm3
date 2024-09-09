#include "Utils.h"
#include "Cell.h"
#include "DataLikelihoodBase.h"
#include "Experiment.h"
#include "NetCDFDataFile.h"
#include "ProbabilityDistributions.h"
#include "SBMLSpecies.h"
#include "Timer.h"
#include "VariabilityDescription.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/random/uniform_01.hpp>
#include <fstream>
#if !PLATFORM_WINDOWS
#include <dlfcn.h>
#endif

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
	, fixed_entry_time(0.0)
	, trailing_simulation_time(0.0)
	, derivative_dll(NULL)
	, simulated_num_cells(0)
	, abs_tol(1e-8)
	, rel_tol(1e-8)
	, requested_duration(EPhaseDuration::None)
	, any_requested_synchronization(false)
	, store_simulation(store_simulation)
	, cell_submit_count(0)
	, cell_done_count(0)
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
	for (size_t i = 0; i < cells.size(); i++) {
		delete cells[i];
	}

	if (derivative_dll) {
#if PLATFORM_WINDOWS
		FreeLibrary(derivative_dll);
#else
		dlclose(derivative_dll);
#endif
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

	if (Cell::use_generated_code) {
		// Generate code for the model/variable combination
		std::string codegen_name = model_filename;
		codegen_name.at(codegen_name.find_last_of('.')) = '_';
		codegen_name += "_" + varset->GetVarsetName() + "_" + Name + "_codegen/";
		if (!GenerateAndCompileSolverCode(codegen_name)) {
			return false;
		}
	}

	// Allocate space for cells
	cells.resize(max_number_of_cells);
	for (size_t i = 0; i < max_number_of_cells; i++) {
		cells[i] = new Cell(&cell_model, this);
	}

	return true;
}

bool Experiment::EvaluateLogProbability(size_t threadix, const VectorReal& values, const VectorReal& transformed_values, Real& logp)
{
	for (size_t i = 0; i < data_likelihoods.size(); i++) {
		data_likelihoods[i]->Reset();
	}
	for (size_t i = 0; i < simulated_trajectories.size(); i++) {
		simulated_cell_parents[i] = std::numeric_limits<size_t>::max();
	}
	if (store_simulation) {
		for (size_t i = 0; i < simulated_trajectories.size(); i++) {
			for (size_t j = 0; j < simulated_trajectories[i].size(); j++) {
				simulated_trajectories[i][j].setConstant(std::numeric_limits<Real>::quiet_NaN());
			}
		}
	}

	active_cells = 0;
	bool result = Simulate(transformed_values);

	if (result) {
		// Report durations
		if (requested_duration != EPhaseDuration::None) {
			for (size_t i = 0; i < active_cells; i++) {
				Real duration = cells[i]->GetDuration(requested_duration);
				for (size_t j = 0; j < data_likelihoods.size(); j++) {
					data_likelihoods[j]->NotifyDuration(i, duration);
				}
			}
		}

		// Synchronize cells if necessary

		// Report data references
		for (int synchronize_i = 0; synchronize_i <= (int)ESynchronizeCellTrajectory::None; synchronize_i++) {
			// Need to reset the timepoint iteration for each different synchronization timepoint
			for (size_t i = 0; i < active_cells; i++) {
				cells[i]->RestartInterpolationIteration();
			}
			for (size_t ti = 0; ti < simulation_timepoints.size(); ti++) {
				const SimulationTimepoints& st = simulation_timepoints[ti];
				if (st.species_ix != std::numeric_limits<size_t>::max() && (int)st.synchronize == synchronize_i) {
					size_t population_size = CountCellsAtTime(st.time, st.synchronize, false);
					size_t mitotic_population_size = CountCellsAtTime(st.time, st.synchronize, true);
					for (size_t i = 0; i < active_cells; i++) {
						Real x = cells[i]->GetInterpolatedSpeciesValue(st.time, st.species_ix, st.synchronize);
						if (x == x) {
							data_likelihoods[st.data_likelihood_ix]->NotifySimulatedValue(st.time_ix, x, st.species_ix, i, population_size, mitotic_population_size, 0, cells[i]->EnteredMitosis(), i >= initial_number_of_cells);
						}
					}
				}
			}
		}

		if (store_simulation) {
			for (size_t i = 0; i < active_cells; i++) {
				cells[i]->RestartInterpolationIteration();
			}
			for (size_t ti = 0; ti < output_trajectories_timepoints.size(); ti++) {
				Real output_t = output_trajectories_timepoints(ti);
				for (size_t i = 0; i < active_cells; i++) {
					for (size_t j = 0; j < cell_model.GetNumCVodeSpecies(); j++) {
						size_t simix = cell_model.GetSimulatedSpeciesFromCVodeSpecies(j);
						simulated_trajectories[i][ti](simix) = cells[i]->GetInterpolatedSpeciesValue(output_t, j, ESynchronizeCellTrajectory::None);
					}
					for (size_t j = 0; j < cell_model.GetNumConstantSpecies(); j++) {
						size_t simix = cell_model.GetSimulatedSpeciesFromConstantSpecies(j);
						simulated_trajectories[i][ti](simix) = cell_model.GetConstantSpecies(j)->GetInitialValue();
					}
					for (size_t j = 0; j < treatment_trajectories.size(); j++) {
						size_t simix = cell_model.GetSimulatedSpeciesByName(cell_model.GetConstantSpeciesName(treatment_trajectories_species_ix[j]));
						simulated_trajectories[i][ti](simix) = treatment_trajectories[j]->GetConcentration(output_t, selected_treatment_trajectory_sample[j]);
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

	Real max_cvode_steps = 0;
	Real min_cvode_step_size = std::numeric_limits<Real>::infinity();

	for (size_t i = 0; i < active_cells; i++) {
		max_cvode_steps = (std::max)(max_cvode_steps, cells[i]->GetCVodeSteps());
		min_cvode_step_size = (std::min)(min_cvode_step_size, cells[i]->GetCVodeMinStepSize());
	}

	cvode_stats.push_back(std::tuple<Real, Real, bool>(max_cvode_steps, min_cvode_step_size, result));

	return true;
}

void Experiment::DumpCVodeStatistics(const std::string& output_folder)
{
	std::string fn = output_folder + std::string("/cvode_statistics_") + Name + std::string(".tsv");
	FILE* f = fopen(fn.c_str(), "w");
	fprintf(f, "cvode_steps\tmin_step_size\tsucceeded\n");
	for (size_t i = 0; i < cvode_stats.size(); i++) {
		fprintf(f, "%g\t%g\t%d\n", std::get<0>(cvode_stats[i]), std::get<1>(cvode_stats[i]), std::get<2>(cvode_stats[i]));
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

	std::string abs_tol_str = xml_node.get<std::string>("<xmlattr>.absolute_tolerance", "");
	std::string rel_tol_str = xml_node.get<std::string>("<xmlattr>.relative_tolerance", "");

	if (abs_tol_str.empty()) {
		abs_tol = 4 * std::numeric_limits<float>::epsilon();
	} else {
		try {
			abs_tol = boost::lexical_cast<Real>(abs_tol_str);
		} catch (const boost::bad_lexical_cast& e) {
			LOGERROR("No conversion possible of absolute tolerance \"%s\" in experiment \"%s\"", abs_tol_str.c_str(), Name.c_str());
			return false;
		}
	}

	if (rel_tol_str.empty()) {
		rel_tol = 4 * std::numeric_limits<float>::epsilon();
	} else {
		try {
			rel_tol = boost::lexical_cast<Real>(rel_tol_str);
		} catch (const boost::bad_lexical_cast& e) {
			LOGERROR("No conversion possible of absolute tolerance \"%s\" in experiment \"%s\"", rel_tol_str.c_str(), Name.c_str());
			return false;
		}
	}

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
			set_init_map[GetCVodeSpeciesByName(name_var.substr(8))] = i;
		}

		if(name_var.substr(0,6).compare("ratio_") == 0){
			// Check that the total variable is also present
			bool also_total_var = false;

			for(int j = 0; j < num_variables; j++){
				std::string name_var_internal = variable_names[j];
				if(name_var_internal.substr(0,6).compare("total_") == 0 && name_var.substr(6).compare(name_var_internal.substr(6)) == 0){
					also_total_var = true;
					std::vector<int> v = {i, j};
					size_t active_species = GetCVodeSpeciesByName("active_" + name_var.substr(6));
					size_t inactive_species = GetCVodeSpeciesByName("inactive_" + name_var.substr(6));

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
				size_t active_species = GetCVodeSpeciesByName("active_" + name_var.substr(6));
				size_t inactive_species = GetCVodeSpeciesByName("inactive_" + name_var.substr(6));

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



#if 0
		// Load treatment trajectories
		auto lapatinib = std::make_unique<TreatmentTrajectory>();
		if (lapatinib->Load(this, &file, "treatment_lapatinib")) {
			treatment_trajectories.push_back(std::move(lapatinib));
			treatment_trajectories_species_ix.push_back(cell_model.GetConstantSpeciesByName("lapatinib"));
			selected_treatment_trajectory_sample.push_back(0);
		} else {
			return false;
		}

		auto trametinib = std::make_unique<TreatmentTrajectory>();
		if (trametinib->Load(this, &file, "treatment_trametinib")) {
			treatment_trajectories.push_back(std::move(trametinib));
			treatment_trajectories_species_ix.push_back(cell_model.GetConstantSpeciesByName("trametinib"));
			selected_treatment_trajectory_sample.push_back(0);
		} else {
			return false;
		}
#endif
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

	Real simulation_begin_time;
	Real simulation_end_time;
	if (!simulation_timepoints.empty()) {
		simulation_begin_time = simulation_timepoints.begin()->time;
		simulation_end_time = simulation_timepoints.rbegin()->time + trailing_simulation_time;
	} else {
		simulation_begin_time = 0;
		simulation_end_time = trailing_simulation_time;
	}

	simulated_cell_parents.resize(max_number_of_cells, std::numeric_limits<size_t>::max());

	if (store_simulation) {
		// Allocate space to store the trajectories
		output_trajectories_timepoints.resize(output_trajectory_num_timepoints);
		for (size_t i = 0; i < output_trajectories_timepoints.size(); i++) {
			output_trajectories_timepoints(i) = simulation_begin_time + (simulation_end_time - simulation_begin_time) * i / (Real)(output_trajectory_num_timepoints - 1);
		}
		simulated_trajectories.resize(max_number_of_cells);
		for (size_t i = 0; i < simulated_trajectories.size(); i++) {
			simulated_trajectories[i].resize(output_trajectory_num_timepoints, VectorReal::Constant(cell_model.GetNumSimulatedSpecies(), std::numeric_limits<Real>::quiet_NaN()));
		}
	}

	return Initialize(xml_node);
}

bool Experiment::Initialize(const boost::property_tree::ptree& xml_node)
{
	std::string entry_time_varname = xml_node.get<std::string>("<xmlattr>.entry_time");
	entry_time_varix = varset->GetVariableIndex(entry_time_varname, false);
	if (entry_time_varix == std::numeric_limits<size_t>::max()) {
		try {
			fixed_entry_time = boost::lexical_cast<Real>(entry_time_varname);
		} catch (const boost::bad_lexical_cast& e) {
			LOGERROR("Could not find variable for mitosis_entry_time \"%s\", and could also not cast it to a constant real value: %s", entry_time_varname.c_str(), e.what());
			return false;
		}
	}

	transformed_variables.setConstant(varset->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());

	if (cell_variabilities.size() > 0) {
		sobol_sequence = std::make_shared< boost::random::sobol >(cell_variabilities.size());
		sobol_sequence_values.resize(initial_number_of_cells * 100);
		boost::random::uniform_01<Real> unif;
		for (int i = 0; i < sobol_sequence_values.size(); i++) {
			sobol_sequence_values[i].resize(cell_variabilities.size());
			for (int j = 0; j < cell_variabilities.size(); j++) {
				sobol_sequence_values[i](j) = unif(*sobol_sequence);
			}
		}
		sobol_sequence_indices.resize(max_number_of_cells, std::numeric_limits<size_t>::max());
	}

	return true;
}

bool Experiment::GenerateAndCompileSolverCode(const std::string& codegen_name)
{
	std::string hash = (boost::format("%x") % std::hash<std::string>{}(codegen_name)).str();
	std::string output_folder = std::string("codegen_") + hash + "/";

	std::string codename_file = output_folder + "codegen_name";
	if (boost::filesystem::exists(codename_file)) {
		std::ifstream infile(codename_file);
		std::string modelfn;
		infile >> modelfn;
		if (modelfn != codegen_name) {
			// Hash mismatch
			LOGERROR("Code generation hashing collision; delete \"%s\" if you wish to regenerate the simulation code", output_folder.c_str());
			return false;
		}
	}

	// Test if the dll is already there.
#if PLATFORM_WINDOWS
#if BUILD_DEBUG
	std::string derivative_dll_fn = output_folder + "Debug/generated_derivatives.dll";
#else
	std::string derivative_dll_fn = output_folder + "Release/generated_derivatives.dll";
#endif
	derivative_dll = LoadLibrary(derivative_dll_fn.c_str());
	if (derivative_dll) {
		derivative = (derivative_fn)GetProcAddress(derivative_dll, "generated_derivative");
		jacobian = (jacobian_fn)GetProcAddress(derivative_dll, "generated_jacobian");
#else
	std::string derivative_dll_fn = output_folder + "libgenerated_derivatives.so";
	derivative_dll = dlopen(derivative_dll_fn.c_str(), RTLD_NOW);
	if (derivative_dll) {
		derivative = (derivative_fn)dlsym(derivative_dll, "generated_derivative");
		jacobian = (jacobian_fn)dlsym(derivative_dll, "generated_jacobian");
#endif
		if (!derivative || !jacobian) {
			LOG("Unable to find generated derivative or jacobian in the dll, recompiling derivative");
			printf("Unable to find generated derivative or jacobian in the dll, recompiling...\n");
		} else {
			LOG("Found DLL with derivative function, reusing it.");
			return true;
		}
	} else {
		LOG("Can't find dll with derivative functions for experiment \"%s\", generating and compiling derivative code.", Name.c_str());
		printf("Can't find dll with derivative functions for experiment \"%s\", generating and compiling derivative code...\n", Name.c_str());
	}

	std::string bcm_path;
	char* bcm_root_env = getenv("BCM3_ROOT");
	if (!bcm_root_env) {
		LOGERROR("BCM3_ROOT environment variable has not been specified!");
		return false;
	} else {
		bcm_path = bcm_root_env;
		bcm3::fix_path(bcm_path);
	}

	std::string code = cell_model.GenerateCode();

	if (!boost::filesystem::is_directory(output_folder)) {
		ASSERT(*output_folder.rbegin() == '/');
		std::string create_path = output_folder.substr(0, output_folder.size() - 1);
		if (!boost::filesystem::create_directories(create_path)) {
			LOGERROR("Unable to make directory for compiling generated code");
			return false;
		}
	}

	std::ofstream f;
	f.open(output_folder + "code.cpp", std::fstream::out);
	if (f.is_open()) {
#if PLATFORM_WINDOWS
		f << "#define WIN32_LEAN_AND_MEAN\n";
		f << "#define NOMINMAX\n";
		f << "#define EIGEN_NO_IO\n";
		f << "#include <windows.h>\n";
#endif
		f << "#include <cmath>\n";
		f << "#include \"sundials/sundials_matrix.h\"\n";
		f << "#include <Eigen/Dense>\n";
		f << "typedef Eigen::VectorXd VectorReal;\n";
		f << "typedef Eigen::MatrixXd MatrixReal;\n";
		f << "#include \"LinearAlgebraSelector.h\"\n";
		f << std::endl;
		for (size_t i = 0; i < cell_model.GetNumCVodeSpecies(); i++) {
			f << "// species[" + std::to_string(i) + "] -- " + cell_model.GetCVodeSpeciesName(i) + " -- " + cell_model.GetCVodeSpecies(i)->GetFullName() + "\n";
		}
		f << std::endl;
		f << "#define EXPORT_PREFIX  extern \"C\"\n";
		f << std::endl;
		f << "inline OdeReal square(OdeReal x) { return x * x; }\n";
		f << std::endl;
		// All of the numerical functions assume that parameters k and n are chosen such that k^n >= epsilon
		// Then x^n + k^n == 0 is due to x^n being 0
		// Also assume that KM's are such that KM^2 >= minimum
		f << "inline OdeReal hill_function(OdeReal x, OdeReal k, OdeReal n)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function(-x, k, n);\n";
		f << "\tif (x == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal xn = pow(x, n);\n";
		f << "\tOdeReal kn = pow(k, n);\n";
		f << "\tOdeReal xnpkn = xn + kn;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn xn / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn1(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function_fixedn1(-x, k);\n";
		f << "\tif (x == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal xnpkn = x + k;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn x / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn2(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function_fixedn2(-x, k);\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tif (x2 == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal xnpkn = x2 + k2;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn x2 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn4(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function_fixedn4(-x, k);\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tif (x4 == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal xnpkn = x4 + k4;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn x4 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn10(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function_fixedn10(-x, k);\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal x8 = x4 * x4;\n";
		f << "\tOdeReal x10 = x2 * x8;\n";
		f << "\tif (x10 == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal k8 = k4 * k4;\n";
		f << "\tOdeReal k10 = k2 * k8;\n";
		f << "\tif (k10 == std::numeric_limits<OdeReal>::infinity()) return 0.0;\n";
		f << "\tOdeReal xnpkn = x10 + k10;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn x10 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn16(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function_fixedn10(-x, k);\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal x8 = x4 * x4;\n";
		f << "\tOdeReal x16 = x8 * x8;\n";
		f << "\tif (x16 == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal k8 = k4 * k4;\n";
		f << "\tOdeReal k16 = k8 * k8;\n";
		f << "\tif (k16 == std::numeric_limits<OdeReal>::infinity()) return 0.0;\n";
		f << "\tOdeReal xnpkn = x16 + k16;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn x16 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn100(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return -hill_function_fixedn100(-x, k);\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal x8 = x4 * x4;\n";
		f << "\tOdeReal x16 = x8 * x8;\n";
		f << "\tOdeReal x32 = x16 * x16;\n";
		f << "\tOdeReal x64 = x32 * x32;\n";
		f << "\tOdeReal x100 = x64 * x32 * x4;\n";
		f << "\tif (x100 == std::numeric_limits<OdeReal>::infinity()) return 1.0;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal k8 = k4 * k4;\n";
		f << "\tOdeReal k16 = k8 * k8;\n";
		f << "\tOdeReal k32 = k16 * k16;\n";
		f << "\tOdeReal k64 = k32 * k32;\n";
		f << "\tOdeReal k100 = k64 * k32 * k4;\n";
		f << "\tif (k100 == std::numeric_limits<OdeReal>::infinity()) return 0.0;\n";
		f << "\tOdeReal xnpkn = x100 + k100;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn x100 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_derivative(OdeReal x, OdeReal k, OdeReal n)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return hill_function_derivative(-x, k, n);\n";
		f << "\tOdeReal xn = pow(x, n);\n";
		f << "\tif (xn == std::numeric_limits<OdeReal>::infinity()) return 0.0;\n";
		f << "\tif (xn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tOdeReal kn = pow(k, n);\n";
		f << "\tif (kn == std::numeric_limits<OdeReal>::infinity()) return 0.0;\n";
		f << "\tOdeReal denom = (x * square(xn + kn));\n";
		f << "\tif (denom < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\treturn kn * n * xn / denom;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_derivative_fixedn1(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x < 0.0) return hill_function_derivative_fixedn1(-x, k);\n";
		f << "\tif (x == std::numeric_limits<OdeReal>::infinity()) return 0.0;\n";
		f << "\tOdeReal xnpkn = x + k;\n";
		f << "\tif (xnpkn == std::numeric_limits<OdeReal>::min()) return 1.0/k;\n";
		f << "\treturn k / (square(xnpkn));\n";
		f << "}\n";
		f << "inline OdeReal michaelis_menten_function(OdeReal kcat, OdeReal KM, OdeReal e, OdeReal s)\n";
		f << "{\n";
		//f << "\tif (s < 0) return -kcat * e * s / (-KM + s);\n";
		f << "\tif (s < 0) return kcat * e * s / KM;\n";
		f << "\treturn kcat * e * s / (KM + s);\n";
		f << "}\n";
		f << "inline OdeReal michaelis_menten_derivative_enzyme(OdeReal kcat, OdeReal KM, OdeReal s)\n";
		f << "{\n";
		//f << "\tif (s < 0) return -kcat * s / (-KM + s);\n";
		f << "\tif (s < 0) return kcat * s / KM;\n";
		f << "\treturn kcat * s / (KM + s);\n";
		f << "}\n";
		f << "inline OdeReal michaelis_menten_derivative_substrate(OdeReal kcat, OdeReal KM, OdeReal e, OdeReal s)\n";
		f << "{\n";
		//f << "\tif (s < 0) return e * kcat * KM / (square(-KM + s));\n";
		f << "\tif (s < 0) return e * kcat / KM;\n";
		f << "\treturn e * kcat * KM / (square(KM + s));\n";
		f << "}\n";
		f << "inline OdeReal safepow(OdeReal x, OdeReal n)\n";
		f << "{\n";
		f << "\tif (x <= 0) {\n";
		f << "\t\treturn 0.0;\n";
		f << "\t} else {\n";
		f << "\t\treturn pow(x, n);\n";
		f << "\t}\n";
		f << "}\n";
		f << std::endl;

		f << code;
		f.close();
	} else {
		LOGERROR("Unable to open output file for generated code");
		return false;
	}

	// Dynamic library definitions file
	f.open(output_folder + "definitions.def", std::fstream::out);
	if (f.is_open()) {
		f << "EXPORTS\n";
		f << "	generated_derivative\n";
		f << "	generated_jacobian\n";
		f.close();
	} else {
		LOGERROR("Unable to open output file for generated code definitions");
		return false;
	}

	f.open(output_folder + "CMakeLists.txt", std::fstream::out);
	if (f.is_open()) {
		f << "cmake_minimum_required(VERSION 3.16)\n";
		f << "project(generated_derivatives CXX)\n";
		f << "\n";
		f << "if(CMAKE_HOST_WIN32)\n";
		f << "	set(CMAKE_CXX_FLAGS \"/arch:AVX2\")\n";
		f << "endif(CMAKE_HOST_WIN32)\n";
		f << "if(CMAKE_HOST_UNIX)\n";
		f << "	set(CMAKE_BUILD_TYPE Release)\n";
		f << "	set(CMAKE_CXX_FLAGS \"-O3 -std=c++11 -march=native\")\n";
		f << "endif(CMAKE_HOST_UNIX)\n";
		f << "\n";
		f << "include_directories(" << bcm_path << "dependencies/eigen-3.4-rc1/ " << bcm_path << "dependencies/cvode-5.3.0/include/ " << bcm_path << "src/odecommon " << bcm_path << "src/utils)\n";
		f << "add_library(generated_derivatives MODULE code.cpp definitions.def)\n";
		f.close();
	} else {
		LOGERROR("Unable to open output file for generated CMake file");
		return false;
	}

	// Compile
	LOG("Compiling generated code...");
	try {
#if PLATFORM_WINDOWS
		// HACKY - assume only one of these is installed and that it's the same one used to compile BCM...
		// Look for vcvars64.bat in the common places
		std::string vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
		FILE *file = fopen(vcvarsfn.c_str(), "r");
		if (file) {
			fclose(file);
#if BUILD_DEBUG
			int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 15 2017 Win64\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
			int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 15 2017 Win64\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
		} else {
			vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
			FILE *file = fopen(vcvarsfn.c_str(), "r");
			if (file) {
				fclose(file);
#if BUILD_DEBUG
				int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 16 2019\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
				int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 16 2019\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
			} else {
				vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
				FILE* file = fopen(vcvarsfn.c_str(), "r");
				if (file) {
					fclose(file);
#if BUILD_DEBUG
					int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
					int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
				} else {
					vcvarsfn = "C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
					FILE* file = fopen(vcvarsfn.c_str(), "r");
					if (file) {
						fclose(file);
#if BUILD_DEBUG
						int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
						int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
					} else {
						LOGERROR("Not sure where to find vcvars64.bat for compiling the code...");
						return false;
					}
				}
			}
		}
#else
		int result = system(("cd " + output_folder + " ; cmake . ; make").c_str());
#endif
	} catch (const std::exception &e) {
		LOGERROR("Process error: %s", e.what());
		return false;
	}

	// Load the compiled dynamic library
#if PLATFORM_WINDOWS
	derivative_dll = LoadLibrary(derivative_dll_fn.c_str());
	if (!derivative_dll) {
		LOGERROR("LoadLibrary failed: %u", GetLastError());
		return false;
	}
	derivative = (derivative_fn)GetProcAddress(derivative_dll, "generated_derivative");
	jacobian = (jacobian_fn)GetProcAddress(derivative_dll, "generated_jacobian");
#else
	derivative_dll = dlopen(derivative_dll_fn.c_str(), RTLD_NOW);
	if (!derivative_dll) {
		LOGERROR("dlopen failed");
		return false;
	}
	derivative = (derivative_fn)dlsym(derivative_dll, "generated_derivative");
	jacobian = (jacobian_fn)dlsym(derivative_dll, "generated_jacobian");
#endif

	if (!derivative || !jacobian) {
		LOGERROR("Unable to find generated derivative or jacobian in the dll");
		return false;
	}

	f.open(codename_file, std::fstream::out);
	if (f.is_open()) {
		f << codegen_name;
		f.close();
	} else {
		LOGERROR("Unable to open codegen name file");
		return false;
	}

	return true;
}

bool Experiment::Simulate(const VectorReal& transformed_values)
{
	for (size_t i = 0; i < varset->GetNumVariables(); i++) {
		transformed_variables[i] = transformed_values[i];
	}
	for (std::map<size_t, size_t>::iterator espi = experiment_specific_parameter_map.begin(); espi != experiment_specific_parameter_map.end(); ++espi) {
		transformed_variables[espi->first] = transformed_variables[espi->second];
	}

	// TODO - incorporate in sobol? Or switch to ML?
	for (size_t i = 0; i < treatment_trajectories.size(); i++) {
		selected_treatment_trajectory_sample[i] = rng.GetUnsignedInt(treatment_trajectories[i]->GetNumSamples()-1);
	}

	Real entry_time;
	if (entry_time_varix != std::numeric_limits<size_t>::max()) {
		entry_time = transformed_variables[entry_time_varix];
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
	Real minimum_start_time = std::numeric_limits<Real>::infinity();

	if (initial_number_of_cells > 1) {
		for (size_t i = 0; i < initial_number_of_cells; i++) {
			Real cell_start_time = -entry_time;
			AddNewCell(cell_start_time, NULL, transformed_values, true, -1);
			minimum_start_time = std::min(cell_start_time, minimum_start_time);
		}
	} else {
		AddNewCell(-entry_time, NULL, transformed_values, false, -1);
		minimum_start_time = -entry_time;
	}

	if (minimum_start_time < -7.0 * 24.0 * 60.0 * 60.0) {
		//printf("Minimum start time too long: %g", minimum_start_time);
		return false;
	}

	for (size_t i = 0; i < cells.size(); i++) {
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

	if (!AuxEvaluationThreads.empty()) {
		cells_to_process_lock.lock();
		aux_target_time = target_time;
		for (size_t i = 0; i < active_cells; i++) {
			cells_to_process.push(i);
		}
		cells_to_process_lock.unlock();

		{
			std::unique_lock<std::mutex> lock(all_done_mutex);
			cell_submit_count = active_cells;
			cell_done_count = 0;
		}

		StartAuxThreads();
		result = WaitAuxThreads();
	} else {
		for (size_t i = 0; i < active_cells; i++) {
			result &= SimulateCell(i, target_time, 0);
			if (!result) {
				break;
			}
		}
	}

	return result;
}

bool Experiment::SimulateCell(size_t i, Real target_time, size_t eval_thread)
{
	bool die = false, divide = false;
	Real achieved_time = 0.0;
	if (!cells[i]->Simulate(target_time, die, divide, achieved_time)) {
		return false;
	}

	if (divide_cells && divide && achieved_time < target_time) {
		size_t cell1, cell2;
		// This is the only place where the number of active cells is increased; the readers don't need to lock it
		{
			std::lock_guard<std::mutex> lock(cell_vector_mutex);

			// Create two new cells
			cell1 = AddNewCell(achieved_time, cells[i], transformed_variables, false, 0);
			cell2 = AddNewCell(achieved_time, cells[i], transformed_variables, false, 1);
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

size_t Experiment::AddNewCell(Real time, Cell* parent, const VectorReal& transformed_values, bool entry_time_variable, int child_ix)
{
	if (active_cells == max_number_of_cells) {
		// TODO - what to do? Fail the simulation or just continue?
		// We'll fail the simulation for now; as it probably indicates the apparent growth rate is too big, at least if the max simulation size was chosen appropriately.
		//LOG("Max number of cells reached\n");
		return std::numeric_limits<size_t>::max();
	}

	// Activate a new cell from the pre-allocated pool
	bool result = true;
	size_t new_cell_ix = active_cells;
	Cell* cell = cells[new_cell_ix];
	active_cells++;

	// Initialize
	cell->SetDerivativeFunctions(derivative, jacobian);

	int sobol_sequence_ix = 0;
	if (parent) {
		result &= cell->SetInitialConditionsFromOtherCell(parent);

		size_t parent_ix = std::find(cells.begin(), cells.end(), parent) - cells.begin();
		simulated_cell_parents[new_cell_ix] = parent_ix;

		if (!sobol_sequence_values.empty()) {
			// This logic assumes max 2 children per cell
			int generation = 0;
			size_t grandparent_ix = simulated_cell_parents[parent_ix];
			while (1) {
				if (grandparent_ix == std::numeric_limits<size_t>::max()) {
					break;
				} else {
					grandparent_ix = simulated_cell_parents[grandparent_ix];
					generation++;
				}
			}
			sobol_sequence_ix = sobol_sequence_indices[parent_ix] * 2 + child_ix;
			while (generation > 0) {
				sobol_sequence_ix += initial_number_of_cells * (1 << generation);
				generation--;
			}
			if (sobol_sequence_ix >= sobol_sequence_values.size()) {
				// Way too many generations, fail the simulation
				return std::numeric_limits<size_t>::max();
			}
		}
	} else {
		result &= cell->SetInitialConditionsFromModel(set_species_map, set_init_map, ratio_active_map, ratio_inactive_map, ratio_total_active, ratio_total_inactive,transformed_values, time);
		if (!sobol_sequence_values.empty()) {
			sobol_sequence_ix = new_cell_ix;
		}
	}
	if (!sobol_sequence_values.empty()) {
		sobol_sequence_indices[new_cell_ix] = sobol_sequence_ix;
	}
	result &= cell->Initialize(time, transformed_values, sobol_sequence_values.empty() ? nullptr : &sobol_sequence_values[sobol_sequence_ix], entry_time_variable, any_requested_synchronization, abs_tol, rel_tol);

	if (!result) {
		return std::numeric_limits<size_t>::max();
	} else {
		return new_cell_ix;
	}
}

size_t Experiment::CountCellsAtTime(Real time, ESynchronizeCellTrajectory synchronize, bool count_only_mitotic)
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

				e.result &= SimulateCell(cell_ix, aux_target_time, threadIndex);
				cells_processed++;

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
