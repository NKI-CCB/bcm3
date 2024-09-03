#include "Utils.h"
#include "SBMLAssignmentRule.h"
#include "SBMLModel.h"
#include "SBMLModelParameters.h"
#include "SBMLReaction.h"
#include "SBMLSpecies.h"

SBMLModel::SBMLModel(size_t numthreads)
	: document(nullptr)
	, model(nullptr)
{
#if BCM3_SBML_INCLUDE_SOLVERS
	Solvers.resize(numthreads);
#endif
}

SBMLModel::~SBMLModel()
{
	delete document;
}

bool SBMLModel::LoadSBML(const std::string& fn)
{
	document = readSBML(fn.c_str());
	if (!document) {
		LOGERROR("Errors reading SBML file %s: NULL returned", fn.c_str());
		return false;
	}
	if (document->getNumErrors(LIBSBML_SEV_ERROR) > 0) {
		LOGERROR("Errors reading SBML file %s: %d", fn.c_str(), document->getNumErrors(LIBSBML_SEV_ERROR));
		for (unsigned int i = 0; i < document->getNumErrors(LIBSBML_SEV_ERROR); i++) {
			LOGERROR(" Error %d: %s", i+1, document->getErrorWithSeverity(i, LIBSBML_SEV_ERROR)->getMessage().c_str());
		}
		delete document;
		document = nullptr;
		return false;
	}
	if (document->getNumErrors(LIBSBML_SEV_FATAL) > 0) {
		LOGERROR("Fatal errors reading SBML file %s: %d", fn.c_str(), document->getNumErrors(LIBSBML_SEV_FATAL));
		for (unsigned int i = 0; i < document->getNumErrors(LIBSBML_SEV_FATAL); i++) {
			LOGERROR(" Error %d: %s", i + 1, document->getErrorWithSeverity(i, LIBSBML_SEV_FATAL)->getMessage().c_str());
		}
		delete document;
		document = nullptr;
		return false;
	}
	if (document->getNumErrors(LIBSBML_SEV_WARNING) > 0) {
		LOGWARNING("Warnings reading SBML file %s: %d", fn.c_str(), document->getNumErrors(LIBSBML_SEV_WARNING));
		for (unsigned int i = 0; i < document->getNumErrors(LIBSBML_SEV_WARNING); i++) {
			LOGWARNING(" Warning %d: %s", i + 1, document->getErrorWithSeverity(i, LIBSBML_SEV_WARNING)->getMessage().c_str());
		}
	}
	if (document->getNumErrors(LIBSBML_SEV_INFO) > 0) {
		LOG("Information while reading SBML file %s: %d", fn.c_str(), document->getNumErrors(LIBSBML_SEV_INFO));
		for (unsigned int i = 0; i < document->getNumErrors(LIBSBML_SEV_INFO); i++) {
			LOG(" Information %d: %s", i + 1, document->getErrorWithSeverity(i, LIBSBML_SEV_INFO)->getMessage().c_str());
		}
	}

	model = document->getModel();
	if (model == nullptr) {
		LOGERROR("Unable to load SBML model");
		delete document;
		document = nullptr;
		return false;
	}

	for (unsigned int i = 0; i < model->getNumSpecies(); i++) {
		std::unique_ptr<SBMLSpecies> sp = std::make_unique<SBMLSpecies>();
		if (!sp->Initialize(model->getSpecies(i), model->getAnnotation())) {
			return false;
		}
		if (Species.find(sp->GetId()) != Species.end()) {
			LOGERROR("Duplicate species id %s", sp->GetId().c_str());
			return false;
		}
		Species[sp->GetId()] = std::move(sp);
	}

	for (unsigned int i = 0; i < model->getNumReactions(); i++) {
		std::unique_ptr<SBMLReaction> re = std::make_unique<SBMLReaction>();
		if (!re->Initialize(model->getReaction(i), this)) {
			return false;
		}
		if (Reactions.find(re->GetId()) != Reactions.end()) {
			LOGERROR("Duplicate reaction id %s", re->GetId().c_str());
			return false;
		}
		Reactions[re->GetId()] = std::move(re);
	}

	// Retrieve the species to simulate
	std::set<std::string> species_without_reactions;
	for (auto si = Species.begin(); si != Species.end(); ++si) {
		if (si->second->GetType() != SBMLSpecies::Sink) {
			SimulatedSpecies.push_back(si->first);
			species_without_reactions.insert(si->first);
		}
	}
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		for (auto rti = ri->second->GetReactants().begin(); rti != ri->second->GetReactants().end(); ++rti) {
			auto si = species_without_reactions.find(*rti);
			if (si != species_without_reactions.end()) {
				species_without_reactions.erase(si);
			}
		}
		for (auto pdi = ri->second->GetProducts().begin(); pdi != ri->second->GetProducts().end(); ++pdi) {
			auto si = species_without_reactions.find(*pdi);
			if (si != species_without_reactions.end()) {
				species_without_reactions.erase(si);
			}
		}
	}
	
	// Do not include species which are neither a reactant nor a product in the CVode calculations
	CVodeSpecies = SimulatedSpecies;
	for (auto ssi = CVodeSpecies.begin(); ssi != CVodeSpecies.end(); ) {
		auto si = species_without_reactions.find(*ssi);
		if (si != species_without_reactions.end()) {
			ConstantSpecies.push_back(*ssi);
			ssi = CVodeSpecies.erase(ssi);
		} else {
			++ssi;
		}
	}
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		if (!ri->second->SetSpeciesVector(CVodeSpecies, ConstantSpecies)) {
			return false;
		}
	}
	cvode_to_simulated_map.resize(CVodeSpecies.size());
	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		auto it = std::find(SimulatedSpecies.begin(), SimulatedSpecies.end(), CVodeSpecies[i]);
		ASSERT(it != SimulatedSpecies.end());
		cvode_to_simulated_map[i] = it - SimulatedSpecies.begin();
	}
	constant_to_simulated_map.resize(ConstantSpecies.size());
	for (size_t i = 0; i < ConstantSpecies.size(); i++) {
		auto it = std::find(SimulatedSpecies.begin(), SimulatedSpecies.end(), ConstantSpecies[i]);
		ASSERT(it != SimulatedSpecies.end());
		constant_to_simulated_map[i] = it - SimulatedSpecies.begin();
	}

	// Initialize assignment rules
	SpeciesIsAssigned.resize(SimulatedSpecies.size(), false);
	for (unsigned int i = 0; i < model->getNumRules(); i++) {
		const Rule* rule = model->getRule(i);
		if (rule->isAssignment()) {
			std::unique_ptr<SBMLAssignmentRule> newrule = std::make_unique<SBMLAssignmentRule>();
			if (!newrule->Initialize(rule, this)) {
				return false;
			}
			newrule->SetSpeciesVector(CVodeSpecies);
			AssignmentRules.push_back(std::move(newrule));

			const std::string& target = newrule->GetVariable();
			bool found = false;
			for (size_t ix = 0; ix < SimulatedSpecies.size(); ix++) {
				if (SimulatedSpecies[ix] == target) {
					AssignmentRulesTargetIx.push_back(ix);
					SpeciesIsAssigned[ix] = true;
					found = true;
					break;
				}
			}

			if (!found) {
				LOGERROR("Target variable \"%s\" for assignment rule %d not found!", target.c_str(), i);
				return false;
			}
		}
	}


	return true;
}

bool SBMLModel::SetVariableSet(std::shared_ptr<const bcm3::VariableSet> variables)
{
	this->variables = variables;

#if BCM3_SBML_INCLUDE_SOLVERS
	for (size_t i = 0; i < Solvers.size(); i++) {
		if (!Solvers[i].parameters->SetVariableSet(variables)) {
			return false;
		}
	}
#endif

	return true;
}

bool SBMLModel::AddNonSampledParameters(const std::vector<std::string>& non_sampled_parameters)
{
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		if (!ri->second->SetNonSampledParameters(non_sampled_parameters)) {
			return false;
		}
	}
	for (auto ri = AssignmentRules.begin(); ri < AssignmentRules.end(); ++ri) {
		if (!(*ri)->SetNonSampledParameters(non_sampled_parameters)) {
			return false;
		}
	}
	return true;
}

bool SBMLModel::AddFixedParameter(const std::string& parameter, Real value)
{
	bool found = false;
	for (unsigned int i = 0; i < model->getNumParameters(); i++) {
		if (model->getParameter(i)->getId() == parameter) {
			found = true;
		}
	}
	if (!found) {
		LOGERROR("Fixed parameter value requested for \"%s\", but there is no parameter with that ID in the SBML model");
		return false;
	}

	fixed_parameter_values[parameter] = value;
	return true;
}

bool SBMLModel::Initialize()
{
	if (!variables) {
		LOGERROR("VariableSet has not been set");
		return false;
	}

	// Map variables to parameters
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		if (!ri->second->SetVariableSet(variables.get())) {
			return false;
		}
		if (!ri->second->PostInitialize(model, fixed_parameter_values)) {
			return false;
		}
	}
	for (size_t ri = 0; ri < AssignmentRules.size(); ri++) {
		if (!AssignmentRules[ri]->SetVariableSet(variables.get())) {
			return false;
		}
		if (!AssignmentRules[ri]->PostInitialize(model, fixed_parameter_values)) {
			return false;
		}
	}

#if BCM3_SBML_INCLUDE_SOLVERS
	// Initialize CVODE solvers
	for (size_t i = 0; i < Solvers.size(); i++) {
		ODESolverCVODE::TDeriviativeFunction derivative = boost::bind(&SBMLModel::CalculateDerivative, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4);

		Solvers[i].parameters = std::make_unique<SBMLModelParameters>(this);
		Solvers[i].solver = std::make_unique<ODESolverCVODE>();
		//Solvers[i].linear_solver = new LinearSolverMKL(this);
		//Solvers[i].solver->SetLinearSolver(Solvers[i].linear_solver);
		Solvers[i].solver->SetDerivativeFunction(derivative);
#if USE_CODE
		ODESolverCVODE::TJacobianFunction jacobian = boost::bind(&SBMLModel::CalculateJacobian, this, boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4, boost::placeholders::_5);
		Solvers[i].solver->SetJacobianFunction(jacobian);
#endif
		if (!Solvers[i].solver->Initialize(CVodeSpecies.size(), Solvers[i].solver.get())) {
			return false;
		}
		Solvers[i].y.resize(CVodeSpecies.size());
		Solvers[i].transformed_variables = VectorReal::Zero(variables->GetNumVariables());
		Solvers[i].transformed_variables_use_in_ode = OdeVectorReal::Zero(variables->GetNumVariables());

		OdeVectorReal tolerance = OdeVectorReal::Zero(CVodeSpecies.size());
		for (size_t si = 0; si < CVodeSpecies.size(); si++) {
			const SBMLSpecies* species = GetSpecies(CVodeSpecies[si]);
			if (species->GetType() == SBMLSpecies::Transcript) {
				tolerance(si) = (OdeReal)1e-10;
			} else {
				tolerance(si) = (OdeReal)1e-6;
			}
		}
		Solvers[i].solver->SetTolerance((OdeReal)1e-5, tolerance);
	}
#endif

	return true;
}

std::string SBMLModel::GenerateCode()
{
	std::string code;

	code += "EXPORT_PREFIX void generated_derivative(OdeReal* out, const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)\n";
	code += "{\n";

	code += "\tOdeReal ratelaws[" + std::to_string(Reactions.size()) + "];\n";
	size_t i = 0;
	std::map<SBMLReaction*, size_t> ratelaw_ix;
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		code += "\tratelaws[" + std::to_string(i) + "] = ";

		std::string eqn = ri->second->GetRateEquation(i);
		if (eqn.empty()) {
			if (ODE_SINGLE_PRECISION) {
				code += "0.0f;\n";
			} else {
				code += "0.0;\n";
			}
		} else {
			code += eqn + ";\n";
		}

		ratelaw_ix[ri->second.get()] = i;
		i++;
	}

	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		code += "\tout[" + std::to_string((uint64)i) + "] = ";

		std::string eqn;
		for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
			const std::vector<std::string>& products = ri->second->GetProducts();
			for (size_t pi = 0; pi < products.size(); pi++) {
				if (products[pi] == CVodeSpecies[i]) {
					Real stoichiometry = ri->second->GetProductStoichiometry(pi);
					if (stoichiometry == 1.0) {
						eqn += "+ratelaws[" + std::to_string(ratelaw_ix[ri->second.get()]) + "]";
					} else if (stoichiometry != 0.0) {
						if (ODE_SINGLE_PRECISION) {
							eqn += "+" + std::to_string(stoichiometry) + "f*ratelaws[" + std::to_string(ratelaw_ix[ri->second.get()]) + "]";
						} else {
							eqn += "+" + std::to_string(stoichiometry) + "*ratelaws[" + std::to_string(ratelaw_ix[ri->second.get()]) + "]";
						}
					}
				}
			}
			const std::vector<std::string>& reactants = ri->second->GetReactants();
			for (size_t pi = 0; pi < reactants.size(); pi++) {
				if (reactants[pi] == CVodeSpecies[i]) {
					Real stoichiometry = ri->second->GetReactantStoichiometry(pi);
					if (stoichiometry == 1.0) {
						eqn += "-ratelaws[" + std::to_string(ratelaw_ix[ri->second.get()]) + "]";
					} else if (stoichiometry != 0.0) {
						if (ODE_SINGLE_PRECISION) {
							eqn += "-" + std::to_string(stoichiometry) + "f*ratelaws[" + std::to_string(ratelaw_ix[ri->second.get()]) + "]";
						} else {
							eqn += "-" + std::to_string(stoichiometry) + "*ratelaws[" + std::to_string(ratelaw_ix[ri->second.get()]) + "]";
						}
					}
				}
			}
		}

		if (eqn.empty()) {
			if (ODE_SINGLE_PRECISION) {
				code += "0.0f;\n";
			} else {
				code += "0.0;\n";
			}
		} else {
			code += eqn + ";\n";
		}
	}

	code += "}\n\n";

	code += "EXPORT_PREFIX void generated_jacobian(SUNMatrix out, const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters)\n";
	code += "{\n";

	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		for (size_t j = 0; j < CVodeSpecies.size(); j++) {
			std::string eqn;
			for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
				ri->second->AddRateDerivativeEquation(i, j, eqn);
			}
			if (!eqn.empty()) {
				code += "\tSM_ELEMENT_D(out, " + std::to_string((uint64)i) + ", " + std::to_string((uint64)j) + ") = " + eqn + ";\n";
			}
		}
		code += "\n";
	}

	code += "}\n";

	return code;
}

std::string SBMLModel::GenerateJacobianCode(size_t i, size_t j)
{
	std::string eqn;
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		ri->second->AddRateDerivativeEquation(j, i, eqn);
	}
	return eqn;
}

bool SBMLModel::CalculateDerivativePublic(OdeReal t, const OdeReal* y, OdeReal* dydt, const OdeReal* constant_species_y, const OdeReal* transformed_variables, const OdeReal* non_sampled_parameters) const
{
	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		dydt[i] = 0.0;
	}
	for (auto ri = Reactions.begin(); ri != Reactions.end(); ++ri) {
		ri->second->AddRates(y, constant_species_y, transformed_variables, non_sampled_parameters, dydt);
	}

	return true;
}

#if BCM3_SBML_INCLUDE_SOLVERS

void SBMLModel::DumpODEStatistics(const std::string& base_filename)
{
	for (size_t i = 0; i < Solvers.size(); i++) {
		std::string fn = base_filename + std::to_string((uint64)i) + std::string(".tsv");
		Solvers[i].solver->DumpStatistics(fn.c_str());
	}
}

void SBMLModel::ResetForcedVariables(size_t threadix)
{
	ParallelSolver& s = Solvers[threadix];
	s.parameters->ResetForcedVariables();
}

bool SBMLModel::SetForcedVariable(const std::string& variable_name, Real value, size_t threadix)
{
	ParallelSolver& s = Solvers[threadix];
	return s.parameters->SetForcedVariable(variable_name, value);
}

bool SBMLModel::Simulate(const Real* variable_values, const Real* initial_conditions, const VectorReal& timepoints, size_t threadix)
{
	ParallelSolver& s = Solvers[threadix];

	// Store variables pointer for use in derivative evaluation
	s.variables = variable_values;
	for (size_t i = 0; i < variables->GetNumVariables(); i++) {
		s.transformed_variables(i) = variables->TransformVariable(i, variable_values[i]);
	}

	// Ensure output is big enough
	if (s.cvode_trajectories.cols() != timepoints.size()) {
		s.cvode_trajectories.setZero(CVodeSpecies.size(), timepoints.size());
	}

	// Force parameters to specific value if mandated by an experimental condition
	s.parameters->UpdateVariableValues(s.transformed_variables.data());

	// Retrieve initial conditions
	if (initial_conditions) {
		for (size_t i = 0; i < CVodeSpecies.size(); i++) {
			s.y[i] = initial_conditions[i];
		}
	} else {
		for (size_t i = 0; i < CVodeSpecies.size(); i++) {
			const SBMLSpecies* sp = GetSpecies(CVodeSpecies[i]);
			std::string param_name = "initial_value_" + sp->GetFullName();
			size_t ix = variables->GetVariableIndex(param_name);
			if (ix == std::numeric_limits<size_t>::max()) {
				s.y[i] = std::numeric_limits<Real>::quiet_NaN();
			} else {
				s.y[i] = s.transformed_variables(ix);
			}
		}
		s.parameters->UpdateInitialValues(s.y.data());
	}
	s.transformed_variables_use_in_ode = s.transformed_variables.cast<OdeReal>();
	CalculateAssignments(s.y.data(), s.transformed_variables_use_in_ode.data(), s.constant_species_y.data(), s.non_sampled_variables.data(), s.y.data());
	
	// Sanity check
	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		if (!(s.y[i] == s.y[i])) {
			LOGERROR("Invalid initial value for species %d/\"%s\"/\"%s\"", i, CVodeSpecies[i].c_str(), GetSimulatedSpeciesFullName(i).c_str());
			LOGERROR("Has an initial value been specified in the prior or as an experimental condition?");
			return false;
		}
	}

	OdeVectorReal timepoints_use_in_solver = timepoints.cast<OdeReal>();
	bool result = s.solver->Simulate(&s.y[0], timepoints_use_in_solver, s.cvode_trajectories);
	if (!result) {
#if 0
		LOGERROR("Model simulation  failed with variable values:");
		for (size_t vi = 0; vi < variables->GetNumVariables(); vi++) {
			LOGERROR("  %s - %g", variables->GetVariableName(vi).c_str(), variable_values[vi]);
		}
#endif
		return false;
	}

	// Fill the final trajectories
	if (s.trajectories.cols() != timepoints.size()) {
		s.trajectories.setZero(SimulatedSpecies.size(), timepoints.size());
	}
	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		size_t ix = GetSimulatedSpeciesByName(CVodeSpecies[i]);
		s.trajectories.row(ix) = s.cvode_trajectories.row(i).cast<Real>();
	}
	for (size_t i = 0; i < ConstantSpecies.size(); i++) {
		// TODO
	}

#if TODO
	// Calculate assignments - NOTE that these assignments are not in effect during the simulation..
	for (size_t ti = 0; ti < timepoints.size(); ti++) {
		VectorReal timepoint = s.trajectories.col(ti);
		CalculateAssignments(timepoint.data(), s.transformed_variables.data(), timepoint.data());
		s.trajectories.col(ti) = timepoint;
	}
#endif

	return true;
}

bool SBMLModel::GetOutput(size_t species_ix, VectorReal& output, size_t threadix)
{
	if (species_ix >= SimulatedSpecies.size()) {
		LOGERROR("Invalid species index");
		return false;
	}

	ParallelSolver& s = Solvers[threadix];
	output = s.trajectories.row(species_ix);
	return true;
}

bool SBMLModel::GetOutput(const std::string& species_full_name, VectorReal& output, size_t threadix)
{
	ParallelSolver& s = Solvers[threadix];

	for (size_t i = 0; i < SimulatedSpecies.size(); i++) {
		if (Species[SimulatedSpecies[i]]->GetFullName() == species_full_name) {
			output = s.trajectories.row(i);
			return true;
		}
	}

	LOGERROR("Could not find species with full name %s", species_full_name.c_str());
	return false;
}

#endif

void SBMLModel::GetParameters(std::vector<std::string>& parameters) const
{
	parameters.resize(model->getNumParameters());
	for (unsigned int i = 0; i < model->getNumParameters(); i++) {
		parameters[i] = model->getParameter(i)->getId();
	}
}

size_t SBMLModel::GetNumSimulatedSpecies() const
{
	return SimulatedSpecies.size();
}

size_t SBMLModel::GetNumCVodeSpecies() const
{
	return CVodeSpecies.size();
}

size_t SBMLModel::GetNumConstantSpecies() const
{
	return ConstantSpecies.size();
}

const SBMLSpecies* SBMLModel::GetSimulatedSpecies(size_t i) const
{
	return Species.at(SimulatedSpecies[i]).get();
}

const SBMLSpecies* SBMLModel::GetCVodeSpecies(size_t i) const
{
	return Species.at(CVodeSpecies[i]).get();
}

const SBMLSpecies* SBMLModel::GetConstantSpecies(size_t i) const
{
	return Species.at(ConstantSpecies[i]).get();
}

const std::string& SBMLModel::GetSimulatedSpeciesName(size_t i) const
{
	return Species.at(SimulatedSpecies[i])->GetName();
}

const std::string& SBMLModel::GetCVodeSpeciesName(size_t i) const
{
	return Species.at(CVodeSpecies[i])->GetName();
}

const std::string& SBMLModel::GetConstantSpeciesName(size_t i) const
{
	return Species.at(ConstantSpecies[i])->GetName();
}

std::string SBMLModel::GetSimulatedSpeciesFullName(size_t i) const
{
	auto iter = Species.find(SimulatedSpecies[i]);
	return iter->second->GetFullName();
}

size_t SBMLModel::GetSimulatedSpeciesByName(const std::string& name, bool log_error) const
{
	for (size_t i = 0; i < SimulatedSpecies.size(); i++) {
		if (Species.at(SimulatedSpecies[i])->GetName() == name) {
			return i;
		}
	}
	if (log_error) {
		LOGERROR("Could not find species \"%s\"", name.c_str());
	}
	return std::numeric_limits<size_t>::max();
}

size_t SBMLModel::GetCVodeSpeciesByName(const std::string& name, bool log_error) const
{
	for (size_t i = 0; i < CVodeSpecies.size(); i++) {
		if (Species.at(CVodeSpecies[i])->GetName() == name) {
			return i;
		}
	}
	if (log_error) {
		LOGERROR("Could not find species \"%s\"", name.c_str());
	}
	return std::numeric_limits<size_t>::max();
}

size_t SBMLModel::GetConstantSpeciesByName(const std::string& name, bool log_error) const
{
	for (size_t i = 0; i < ConstantSpecies.size(); i++) {
		if (Species.at(ConstantSpecies[i])->GetName() == name) {
			return i;
		}
	}
	if (log_error) {
		LOGERROR("Could not find species \"%s\"", name.c_str());
	}
	return std::numeric_limits<size_t>::max();
}

void SBMLModel::GetSimulatedSpeciesFullNames(std::vector<std::string>& species_names) const
{
	species_names.clear();
	for (size_t i = 0; i < SimulatedSpecies.size(); i++) {
		auto iter = Species.find(SimulatedSpecies[i]);
		species_names.push_back(iter->second->GetFullName());
	}
}

const SBMLSpecies* SBMLModel::GetSpecies(const std::string& id) const
{
	auto iter = Species.find(id);
	if (iter != Species.end()) {
		return iter->second.get();
	} else {
		return NULL;
	}
}

const SBMLSpecies* SBMLModel::GetSpeciesFromFullName(const std::string& id) const
{
	for (auto si = Species.begin(); si != Species.end(); ++si) {
		if (si->second->GetFullName() == id) {
			return si->second.get();
		}
	}
	return NULL;
}

#if BCM3_SBML_INCLUDE_SOLVERS

bool SBMLModel::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user) const
{
	// Find which parallel solver was used
	size_t threadix;
	for (size_t i = 0; i < Solvers.size(); i++) {
		if (Solvers[i].solver.get() == user) {
			threadix = i;
			break;
		}
	}

	// Update constant species

#if USE_CODE
	derivative(dydt, y, Solvers[threadix].transformed_variables_use_in_ode.data());
	return true;
#else
	return CalculateDerivativePublic(t, y, dydt, Solvers[threadix].constant_species_y.data(), Solvers[threadix].transformed_variables_use_in_ode.data(), Solvers[threadix].non_sampled_variables.data());
#endif
}

bool SBMLModel::CalculateJacobian(OdeReal t, const OdeReal* y, const OdeReal* dydt, OdeReal** jac, void* user) const
{
	// Find which parallel solver was used
	size_t threadix;
	for (size_t i = 0; i < Solvers.size(); i++) {
		if (Solvers[i].solver.get() == user) {
			threadix = i;
			break;
		}
	}
	
#if USE_CODE
	jacobian(jac, y, Solvers[threadix].transformed_variables_use_in_ode.data());
#endif

	return true;
}

void SBMLModel::CalculateAssignments(const OdeReal* y, const OdeReal* variables, const OdeReal* constant_species_y, const OdeReal* non_sampled_parameters, OdeReal* out) const
{
	for (size_t ri = 0; ri < AssignmentRules.size(); ri++) {
		OdeReal value;
		AssignmentRules[ri]->Calculate(y, constant_species_y, variables, non_sampled_parameters, &value);
		out[AssignmentRulesTargetIx[ri]] = value;
	}
}

#endif
