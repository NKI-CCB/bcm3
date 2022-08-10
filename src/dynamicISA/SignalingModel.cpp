#include "Utils.h"
#include "SignalingLink.h"
#include "SignalingModel.h"

#include <sbml/SBMLTypes.h>
using namespace boost::placeholders;

SignalingModel::SolverParams::SolverParams()
	: base_param_values(nullptr)
	, strength_param_values(nullptr)
	, inhib_param_values(nullptr)
	, expression_mixing_param_values(nullptr)
	, fixed_species_values(nullptr)
	, expression_levels(nullptr)
{
}

SignalingModel::SignalingModel()
{
}

SignalingModel::~SignalingModel()
{
	solver.DumpStatistics("cvode_statistics.log");
}

bool SignalingModel::Load(const std::string& sbml_fn)
{
	// Load SBML file
	SBMLDocument* document = readSBML(sbml_fn.c_str());
	if (!document) {
		LOGERROR("Errors reading SBML file %s: NULL return value", sbml_fn.c_str());
		return false;
	}
	if (document->getNumErrors() > 0) {
		LOGERROR("Errors reading SBML file %s: %d", sbml_fn.c_str(), document->getNumErrors());
		for (unsigned int i = 0; i < document->getNumErrors(); i++) {
			LOGERROR(" %d: %s", i, document->getError(i)->getMessage().c_str());
		}
		delete document;
		return false;
	}

	const Model* model = document->getModel();

	for (unsigned int i = 0; i < model->getNumSpecies(); i++) {
		const Species* sp = model->getSpecies(i);
		std::unique_ptr<SignalingMolecule> m = SignalingMolecule::Create(sp);
		if (m) {
			signaling_molecules.push_back(std::move(m));
		} else {
			delete document;
			return false;
		}
	}

	for (unsigned int i = 0; i < model->getNumReactions(); i++) {
		const Reaction* re = model->getReaction(i);
		std::unique_ptr<SignalingLink> l = SignalingLink::Create(re, this);
		if (!l) {
			delete document;
			return false;
		}

		bool source_exists = false;
		bool target_exists = false;
		for (auto smi = signaling_molecules.begin(); smi != signaling_molecules.end(); smi++) {
			if ((*smi)->GetName() == l->GetFrom()) {
				source_exists = true;
			}
			if ((*smi)->GetName() == l->GetTo()) {
				target_exists = true;
			}
		}

		if (!source_exists || !target_exists) {
			LOGERROR("Signaling link from \"%s\" to \"%s\" references unknown species", l->GetFrom().c_str(), l->GetTo().c_str());
			delete document;
			return false;
		}

		signaling_links.push_back(std::move(l));
	}

	std::set<std::string> simulated_species_names;
	for (size_t i = 0; i < signaling_links.size(); i++) {
		simulated_species_names.insert(signaling_links[i]->GetTo());
	}

	simulated_species.resize(simulated_species_names.size());
	std::vector< std::string > simulated_species_names_vector(simulated_species_names.size());
	size_t i = 0;
	for (auto ssni = simulated_species_names.begin(); ssni != simulated_species_names.end(); ssni++, i++) {
		size_t ix = 0;
		for (auto smi = signaling_molecules.begin(); smi != signaling_molecules.end(); smi++, ix++) {
			if ((*smi)->GetName() == *ssni) {
				break;
			}
		}
		simulated_species[i] = ix;
		simulated_species_names_vector[i] = *ssni;
	}

	for (auto sli = signaling_links.begin(); sli != signaling_links.end(); ++sli) {
		(*sli)->NotifySimulatedSpecies(simulated_species_names_vector);
	}

	CVODESolver::TDeriviativeFunction derivative = boost::bind(&SignalingModel::CalculateDerivative, this, _1, _2, _3, _4);
	solver.Initialize(simulated_species.size(), this);
	solver.SetDerivativeFunction(derivative);
	solver.SetTolerance(1e-6, 1e-6);

	return true;
}

bool SignalingModel::Calculate(const VectorReal& base_param_values,
							   const VectorReal& decay_param_values,
							   const VectorReal& strength_param_values,
							   const VectorReal& inhib_param_values,
							   const VectorReal& expression_mixing_param_values,
							   const VectorReal& fixed_species_values,
							   const VectorReal& expression_levels,
							   Real cell_cycle_duration,
							   VectorReal& activities)
{
	solver_params.base_param_values = &base_param_values;
	solver_params.decay_param_values = &decay_param_values;
	solver_params.strength_param_values = &strength_param_values;
	solver_params.inhib_param_values = &inhib_param_values;
	solver_params.expression_mixing_param_values = &expression_mixing_param_values;
	solver_params.expression_levels = &expression_levels;

	if (!SetupCalculations(base_param_values, fixed_species_values)) {
		return false;
	}

	VectorReal timepoints(20);
	for (unsigned int i = 0; i < timepoints.size(); i++) {
		timepoints(i) = i * cell_cycle_duration / ((Real)timepoints.size() - 1.0);
	}
	VectorReal initial_conditions = VectorReal::Zero(simulated_species.size());
	MatrixReal output(simulated_species.size(), timepoints.size());
	if (!solver.Simulate(initial_conditions.data(), timepoints, output)) {
		return false;
	}

	activities = solver_params.constant_species_values;
	for (index_type i = 0; i < simulated_species.size(); i++) {
		index_type ix = simulated_species[i];
		activities(ix) = 0.0;
		for (size_t j = 1; j < timepoints.size(); j++) {
			Real dt = timepoints(j) - timepoints(j-1);
			Real x = 0.5 * (output(i,j) + output(i, j-1));
			activities(ix) += x * dt;
		}
		activities(ix) *= 1.0 / cell_cycle_duration;
	}

	return true;
}

bool SignalingModel::Calculate(const VectorReal& base_param_values,
	const VectorReal& decay_param_values,
	const VectorReal& strength_param_values,
	const VectorReal& inhib_param_values,
	const VectorReal& expression_mixing_param_values,
	const VectorReal& fixed_species_values,
	const VectorReal& expression_levels,
	VectorReal& timepoints,
	std::vector< std::tuple<size_t, Real, Real> >& treatments,
	std::vector<VectorReal>& activities)
{
	solver_params.base_param_values = &base_param_values;
	solver_params.decay_param_values = &decay_param_values;
	solver_params.strength_param_values = &strength_param_values;
	solver_params.inhib_param_values = &inhib_param_values;
	solver_params.expression_mixing_param_values = &expression_mixing_param_values;
	solver_params.expression_levels = &expression_levels;
	solver_params.treatments = &treatments;
	solver_params.current_discontinuity_ix = 0;

	if (!SetupCalculations(base_param_values, fixed_species_values)) {
		return false;
	}
	if (!treatments.empty()) {
		solver.SetDiscontinuity(std::get<1>((*solver_params.treatments)[0]), boost::bind(&SignalingModel::TreatmentCallback, this, _1, _2), nullptr);
	} else {
		solver.ResetDiscontinuity();
	}

	VectorReal initial_conditions = VectorReal::Zero(simulated_species.size());
	MatrixReal output(simulated_species.size(), timepoints.size());
	if (!solver.Simulate(initial_conditions.data(), timepoints, output)) {
		return false;
	}

	for (size_t j = 0; j < timepoints.size(); j++) {
		for (index_type i = 0; i < simulated_species.size(); i++) {
			index_type ix = simulated_species[i];
			activities[j](ix) = output(i,j);
		}
	}

	return true;
}

bool SignalingModel::SetupCalculations(const VectorReal& base_param_values, const VectorReal& fixed_species_values)
{
	solver_params.constant_species_values.setZero(signaling_molecules.size());
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		if (!std::isnan(fixed_species_values(i))) {
			solver_params.constant_species_values(i) = fixed_species_values(i);
		} else {
			solver_params.constant_species_values(i) = base_param_values(i);
		}
	}

	if (!CalculateInhibitions()) {
		return false;
	}
	return true;
}

bool SignalingModel::CalculateInhibitions()
{
	solver_params.inhibitions.setOnes(signaling_molecules.size());
	for (size_t i = 0; i < signaling_links.size(); i++) {
		const SignalingLink* sl = signaling_links[i].get();
		const SignalingMolecule* sm = signaling_molecules[sl->GetFromIx()].get();
		if (sm->type == SignalingMolecule::Drug) {
			ASSERT(sl->FromIsConstant());
			Real inhibition = 1.0 - solver_params.constant_species_values[sl->GetFromIx()] * (*solver_params.strength_param_values)(i);
			solver_params.inhibitions(sl->GetToIx()) *= inhibition;
		}
	}

	return true;
}

bool SignalingModel::CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user)
{
	for (index_type i = 0; i < simulated_species.size(); i++) {
		index_type ix = simulated_species[i];
		dydt[i] = (*solver_params.base_param_values)[ix];
	}

	for (size_t i = 0; i < signaling_links.size(); i++) {
		const SignalingLink* sl = signaling_links[i].get();
		const SignalingMolecule* sm = signaling_molecules[sl->GetFromIx()].get();
		if (sm->type == SignalingMolecule::Drug) {
			continue;
		}

		Real factor;
		if (sl->IsActivating()) {
			factor = 1.0;
		} else {
			factor = -1.0;
		}

		if (sl->FromIsConstant()) {
			dydt[sl->GetToSimulatedIx()] += factor * solver_params.constant_species_values[sl->GetFromIx()] * (*solver_params.strength_param_values)(i);
		} else {
			dydt[sl->GetToSimulatedIx()] += factor * y[sl->GetFromSimulatedIx()] * (*solver_params.strength_param_values)(i);
		}
	}

	for (index_type i = 0; i < simulated_species.size(); i++) {
		index_type ix = simulated_species[i];
		if (dydt[i] > 0) {
			dydt[i] *= (1.0 - y[i]);
		} else {
			dydt[i] *= y[i];
		}
		dydt[i] *= (*solver_params.expression_levels)(ix) * (*solver_params.expression_mixing_param_values)(ix) + (1.0 - (*solver_params.expression_mixing_param_values)(ix));
		dydt[i] *= solver_params.inhibitions(ix);
		dydt[i] -= (*solver_params.decay_param_values)(ix) * y[i];
	}

	return true;
}

Real SignalingModel::TreatmentCallback(OdeReal t, void* user)
{
	auto tr = (*solver_params.treatments)[solver_params.current_discontinuity_ix];
	size_t ix = std::get<0>(tr);
	solver_params.constant_species_values[ix] = std::get<2>(tr);

	CalculateInhibitions();

	solver_params.current_discontinuity_ix++;
	if (solver_params.current_discontinuity_ix < solver_params.treatments->size()) {
		return std::get<1>((*solver_params.treatments)[solver_params.current_discontinuity_ix]);
	} else {
		return std::numeric_limits<Real>::infinity();
	}
}

std::string SignalingModel::GetMoleculeNameById(const std::string& id) const
{
	for (index_type i = 0; i < signaling_molecules.size(); i++) {
		if (signaling_molecules[i]->GetId() == id) {
			return signaling_molecules[i]->GetName();
		}
	}
	LOGERROR("Signaling molecule id \"%s\" not found", id.c_str());
	return std::string();
}

size_t SignalingModel::GetMoleculeIxByName(const std::string& name, bool log_error) const
{
	for (index_type i = 0; i < signaling_molecules.size(); i++) {
		if (signaling_molecules[i]->GetName() == name) {
			return i;
		}
	}
	if (log_error) {
		LOGERROR("Signaling molecule \"%s\" not found", name.c_str());
	}
	return std::numeric_limits<size_t>::max();
}
