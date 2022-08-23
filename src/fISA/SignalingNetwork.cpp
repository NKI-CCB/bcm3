#include "Utils.h"
#include "SignalingNetwork.h"

#include <sbml/SBMLTypes.h>
#include <boost/algorithm/string.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/random/uniform_01.hpp>

static const Real FIXED_K = 9.19024;

inline Real square(Real a) { return a*a; }
inline Real logistic_activation(Real x, Real steepness, Real inflection)
{
	return 1.0 / (1.0 + exp(-steepness * (x - inflection)));
}
inline Real logistic_activation_fixed_bks(Real sum)
{
	if (sum > 3.5) {
		return 1.0;
	} else {
		return 1.0 / (1.0 + exp(-9.19024 * (sum - 0.5)));
	}
}

inline Real logistic_deriv_part(Real sum)
{
	//if (sum < -3.5 || sum > 3.5) {
	//	return 0.0;
	//} else {
		Real g = exp(-FIXED_K * (sum - 0.5));
		return g / square(g + 1.0);
	//}
}

Real SignalingNetwork::max_expression_function(Real expression, size_t signaling_molecule_ix, const Real* values) const
{
	const SignalingMolecule& sm = signaling_molecules[signaling_molecule_ix];
	return expression_function(1.0, expression, sm.expression_mixing_param, values);
}

inline Real SignalingNetwork::expression_function(Real activity, Real expression, SignalingNetwork::parameter_index_type expression_mixing_param, const Real* values) const
{
	if (expression_mixing_param != std::numeric_limits<SignalingNetwork::parameter_index_type>::max()) {
		Real em = GetValue(values, expression_mixing_param);
		return (em * expression + (1.0 - em)) * activity;
	} else {
		return expression * activity;
	}
}

SignalingNetwork::SignalingMolecule::SignalingMolecule()
	: type(InvalidType)
	, drug_inhibition_type(InvalidDrugInhibitionType)
	, base_parameter_ix(std::numeric_limits<parameter_index_type>::max())
	, num_parents(0)
	, input_mixing_param(Mixing_Undefined)
	, expression_mixing_param(std::numeric_limits<parameter_index_type>::max())
	, has_susceptibility_altering_parent(false)
{
	for (int i = 0; i < MAX_PARENTS; i++) {
		drug_susceptibility[i] = std::numeric_limits<parameter_index_type>::max();
	}
}

SignalingNetwork::DrugTransporter::DrugTransporter()
	: num_targets(0)
{
}

SignalingNetwork::SignalingNetwork(size_t numthreads)
	: numthreads(numthreads)
	, VarSet(NULL)
	, activation_limit(ActivationLimit_Invalid)
	, num_strongly_connected_components(0)
{
}

SignalingNetwork::~SignalingNetwork()
{
}

bool SignalingNetwork::Initialize(const std::string& sbml_file, std::shared_ptr<const bcm3::VariableSet> varset, ActivationLimitType activation_limit, size_t multiroot_solves)
{
	if (activation_limit != ActivationLimit_MinMax && activation_limit != ActivationLimit_Logistic && activation_limit != ActivationLimit_LogisticOr) {
		LOGERROR("Invalid activation limit");
		return false;
	}

	VarSet = varset;
	this->activation_limit = activation_limit;
	this->multiroot_solves = multiroot_solves;

	SBMLDocument* document = NULL;
	const Model* model = NULL;

	// Load SBML file
	document = readSBML(sbml_file.c_str());
	if (!document) {
		LOGERROR("Errors reading SBML file %s: NULL return value");
		return false;
	}
	if (document->getNumErrors() > 0) {
		LOGERROR("Errors reading SBML file %s: %d", sbml_file.c_str(), document->getNumErrors());
		for (unsigned int i = 0; i < document->getNumErrors(); i++) {
			LOGERROR(" %d: %s", i, document->getError(i)->getMessage().c_str());
		}
		delete document;
		return false;
	}

	model = document->getModel();

	for (unsigned int i = 0; i < model->getNumSpecies(); i++) {
		const Species* sp = model->getSpecies(i);
		const XMLNode* annotation = sp->getAnnotation();

		SignalingMolecule s;
		s.id = sp->getId();
		s.name = sp->getName();

		for (unsigned int ai = 0; ai < annotation->getNumChildren(); ai++) {
			const XMLNode& child = annotation->getChild(ai);
			if (child.getPrefix() == "celldesigner" && child.getName() == "extension") {
				std::string species_class = child.getChild("speciesIdentity").getChild("class").getChild(0).getCharacters();
				if (species_class.empty()) {
					LOGERROR("Can't find species class for species %s [%s] - SBML file should be a CellDesigner file", sp->getId().c_str(), sp->getName().c_str());
					delete document;
					return false;
				}

				if (species_class == "PROTEIN") {
					s.type = SignalingMolecule::Protein;
				} else if (species_class == "RNA") {
					s.type = SignalingMolecule::mRNA;
				} else if (species_class == "SIMPLE_MOLECULE") {
					s.type = SignalingMolecule::SmallMolecule;
				} else if (species_class == "GENE") {
					s.type = SignalingMolecule::Mutation;
				} else if (species_class == "DRUG") {
					s.type = SignalingMolecule::Drug;
				} else if (species_class == "PHENOTYPE") {
					s.type = SignalingMolecule::Phenotype;
				} else if (species_class == "UNKNOWN") {
					s.type = SignalingMolecule::Unknown;
				} else {
					LOGERROR("Unrecognized species type %s for species %s [%s]", species_class.c_str(), sp->getId().c_str(), sp->getName().c_str());
					delete document;
					return false;
				}
			}
		}

		std::string notes;
		const XMLNode* notes_node = sp->getNotes();
		if (notes_node) {
			const XMLNode& child = notes_node->getChild(0);
			const XMLNode& child1 = child.getChild(1);
			const XMLNode& child2 = child1.getChild(0);
			notes = child2.getCharacters();
			boost::trim(notes);
		}

		if (s.type == SignalingMolecule::Drug) {
			if (notes == "inhibit activity") {
				s.drug_inhibition_type = SignalingMolecule::InhibitActivity;
			} else if (notes == "inhibit activity,alter susceptibility") {
				s.drug_inhibition_type = SignalingMolecule::InhibitActivityAlterSusceptibility;
			} else if (notes == "alter susceptibility") {
				s.drug_inhibition_type = SignalingMolecule::AlterSusceptibility;
			} else if (notes == "inhibit activation") {
				s.drug_inhibition_type = SignalingMolecule::InhibitActivation;
			} else if (notes == "activate") {
				s.drug_inhibition_type = SignalingMolecule::Activate;
			} else if (notes == "") {
				LOGERROR("The drug \"%s\" does not have a note to specify whether the drug inhibits the activity or the activation of its targets. Please add a \"note\" to the species in the SBML file specifying either \"inhibit activity\", \"inhibit activation\" or \"activate\".", s.name.c_str());
				delete document;
				return false;
			} else {
				LOGERROR("The drug \"%s\" has notes, but the notes do not appear to specify whether the drug inhibits the activity or the activation of its targets. Please add a \"note\" to the species in the SBML file specifying either \"inhibit activity\", \"inhibit activation\" or \"activate\".", s.name.c_str());
				delete document;
				return false;
			}
		} else {
			s.drug_inhibition_type = SignalingMolecule::NotADrug;

			if (s.type == SignalingMolecule::Protein) {
				if (notes == "drug_transporter") {
					s.type = SignalingMolecule::DrugTransporter;

					DrugTransporter dt;
					dt.id = s.id;
					dt.name = s.name;
					dt.smix = (signaling_index_type)signaling_molecules.size();
					drug_transporters.push_back(dt);
				}
			} else if (s.type == SignalingMolecule::Mutation) {
				if (notes == "complete_loss") {
					s.type = SignalingMolecule::CompleteLossMutation;
				}
			}
		}

		if (notes == "input_mixing_multiply") {
			s.input_mixing_param = SignalingMolecule::Mixing_Multiply;
		}
		if (notes == "input_mixing_sum") {
			s.input_mixing_param = SignalingMolecule::Mixing_Sum;
		}

		size_t expression_mixing_param = VarSet->GetVariableIndex(std::string("expression_mixing[") + s.name + std::string("]"), false);
		if (expression_mixing_param != std::numeric_limits<size_t>::max()) {
			s.expression_mixing_param = expression_mixing_param;
		}

		signaling_molecules.push_back(s);
	}

	for (unsigned int i = 0; i < model->getNumReactions(); i++) {
		const Reaction* re = model->getReaction(i);

		const XMLNode* annotation = re->getAnnotation();
		bool activating = true;
		for (unsigned int i = 0; i < annotation->getNumChildren(); i++) {
			const XMLNode& child = annotation->getChild(i);
			if (child.getPrefix() == "celldesigner" && child.getName() == "extension") {
				std::string reaction_type_str = child.getChild("reactionType").getChild(0).getCharacters();
				if (reaction_type_str.empty()) {
					LOGERROR("Can't find reaction type for reaction %s - SBML file should be a CellDesigner file", re->getId().c_str());
					delete document;
					return false;
				}
			
				if (reaction_type_str == "POSITIVE_INFLUENCE") {
					activating = true;
				} else if (reaction_type_str == "NEGATIVE_INFLUENCE") {
					activating = false;
				} else {
					LOGERROR("Unrecognized reaction type %s for reaction %s - SBML file should use reduced notation", reaction_type_str.c_str(), re->getId().c_str());
					delete document;
					return false;
				}
			}
		}

		if (re->getNumReactants() != 1) {
			LOGERROR("Can only handle reactions with 1 reactant");
			delete document;
			return false;
		}
		if (re->getNumProducts() != 1) {
			LOGERROR("Can only handle reactions with 1 product");
			delete document;
			return false;
		}

		size_t product_ix = GetSignalingMoleculeIxById(re->getProduct(0)->getSpecies());
		if (product_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Unknown product \"%s\"", re->getProduct(0)->getSpecies().c_str());
			delete document;
			return false;
		}

		if (signaling_molecules[product_ix].type == SignalingMolecule::Drug) {
			// This is a drug transporter reaction
			std::string reactant_id = re->getReactant(0)->getSpecies();
			bool found = true;
			for (size_t dti = 0; dti < drug_transporters.size(); dti++) {
				DrugTransporter& dt = drug_transporters[dti];
				if (dt.id == reactant_id) {
					dt.drug_target[dt.num_targets] = (signaling_index_type)product_ix;
					dt.import[dt.num_targets] = activating ? 1 : 0;
					if (!GetVariableIx(varset.get(), std::string("strength_") + dt.name + std::string("_") + signaling_molecules[product_ix].name, dt.strength_varix[dt.num_targets])) {
						return false;
					}
					dt.num_targets++;
					break;
				}
			}

			if (!found) {
				LOGERROR("Reaction affects a drug, but the reactant has not been marked as drug transporter.");
				return false;
			}
		} else {
			// Normal reaction
			size_t reactant_ix = GetSignalingMoleculeIxById(re->getReactant(0)->getSpecies());
			if (reactant_ix == std::numeric_limits<size_t>::max()) {
				LOGERROR("Unknown reactant \"%s\"", re->getReactant(0)->getSpecies().c_str());
				delete document;
				return false;
			}

			if (signaling_molecules[product_ix].num_parents >= MAX_PARENTS) {
				LOGERROR("Node \"%s\" has too many incoming reactions", signaling_molecules[product_ix].name.c_str());
				return false;
			}
			if (reactant_ix >= std::numeric_limits<signaling_index_type>::max()) {
				LOGERROR("Network too large");
				return false;
			}
			signaling_index_type parent_ix = signaling_molecules[product_ix].num_parents;
			signaling_molecules[product_ix].num_parents++;

			signaling_molecules[product_ix].parents[parent_ix] = (signaling_index_type)reactant_ix;
			signaling_molecules[product_ix].parents_activating[parent_ix] = activating;

			signaling_molecules[product_ix].param1[parent_ix] = std::numeric_limits<parameter_index_type>::max();
			signaling_molecules[product_ix].param2[parent_ix] = std::numeric_limits<parameter_index_type>::max();
			signaling_molecules[product_ix].param3[parent_ix] = std::numeric_limits<parameter_index_type>::max();

			// Drugs are non-linear by default
			if (signaling_molecules[reactant_ix].type == SignalingMolecule::Drug) {
				signaling_molecules[product_ix].param2[parent_ix] = 1;
			}

			// Check for non-linear reactions
			const XMLNode* notes = re->getNotes();
			if (notes) {
				const XMLNode& child = notes->getChild(0);
				const XMLNode& child1 = child.getChild(1);
				const XMLNode& child2 = child1.getChild(0);

				std::string reaction_notes = child2.getCharacters();
				boost::trim(reaction_notes);

				if (reaction_notes == "nonlinear") {
					// Temporarily set to 1 so we know to look it up later.
					signaling_molecules[product_ix].param2[parent_ix] = 1;
				}
				else if (reaction_notes == "linear") {
					signaling_molecules[product_ix].param2[parent_ix] = std::numeric_limits<parameter_index_type>::max();
				}
			}
		}

		// TODO - check for duplicate reactions?
	}

	// Set up variable mapping
	bool result = true;
	signaling_index_type drug_inhibition_counter = 0;
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		SignalingMolecule& sm = signaling_molecules[i];
		
		if (sm.type != SignalingMolecule::Mutation && sm.type != SignalingMolecule::CompleteLossMutation && sm.type != SignalingMolecule::Drug && sm.type != SignalingMolecule::DrugTransporter) {
			if (!GetVariableIx(varset.get(), std::string("base_") + sm.name, sm.base_parameter_ix)) {
				LOGWARNING("No variable for base activation level of %s - assuming 0 base activity unless likelihood specifies the activity.", sm.name.c_str());
			}
		}

		//if (sm.type == SignalingMolecule::Drug && sm.name == drug) {
		//	drug_signaling_ix = i;
		//}
		
		for (signaling_index_type j = 0; j < sm.num_parents; j++) {
			const SignalingMolecule& parent = signaling_molecules[sm.parents[j]];

			// Look up parameters for the activation function (can be linear or nonlinear)
			if (parent.type == SignalingMolecule::Drug) {
				if (parent.drug_inhibition_type == SignalingMolecule::InhibitActivity ||
					parent.drug_inhibition_type == SignalingMolecule::InhibitActivityAlterSusceptibility ||
					parent.drug_inhibition_type == SignalingMolecule::InhibitActivation ||
					parent.drug_inhibition_type == SignalingMolecule::Activate) {
					result &= GetVariableIx(varset.get(), std::string("maxinhib_") + parent.name + std::string("_") + sm.name, sm.param1[j]);
					if (sm.param2[j] == 1) {
						result &= GetVariableIx(varset.get(), std::string("ic50_") + parent.name + std::string("_") + sm.name, sm.param2[j]);
						result &= GetVariableIx(varset.get(), std::string("logsteepness_") + parent.name + std::string("_") + sm.name, sm.param3[j]);
					}
					sm.drug_inhibition_ix[j] = drug_inhibition_counter++;
				}
			} else if (parent.type == SignalingMolecule::CompleteLossMutation) {
				// No parameters.
			} else {
				result &= GetVariableIx(varset.get(), std::string("strength_") + parent.name + std::string("_") + sm.name, sm.param1[j]);
				if (sm.param2[j] == 1) {
					result &= GetVariableIx(varset.get(), std::string("inflection_") + parent.name + std::string("_") + sm.name, sm.param2[j]);
					result &= GetVariableIx(varset.get(), std::string("steepness_") + parent.name + std::string("_") + sm.name, sm.param3[j]);
				}
			}

			// Are there any parent drugs which alter susceptibility?
			if (parent.type == SignalingMolecule::Drug && (parent.drug_inhibition_type == SignalingMolecule::InhibitActivityAlterSusceptibility ||
														   parent.drug_inhibition_type == SignalingMolecule::AlterSusceptibility)) {
				result &= GetVariableIx(varset.get(), parent.name + "_" + sm.name + "_susceptibility", sm.drug_susceptibility[j]);
				if (sm.has_susceptibility_altering_parent) {
					LOGERROR("Error: signaling molecule has more than 1 drug which alters susceptibility - this is currently not implemented");
				}
				sm.has_susceptibility_altering_parent = true;
			}
		}
	}
	if (!result) {
		delete document;
		return false;
	}
	
	// Calculate number of inputs factors
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		SignalingMolecule& sm = signaling_molecules[i];

		size_t activating_input_count = 0;
		size_t inhibiting_input_count = 0; // excl drugs
		
		for (size_t j = 0; j < sm.num_parents; j++) {
			size_t parent_ix = sm.parents[j];
			const SignalingMolecule& parent = signaling_molecules[parent_ix];
			if (parent.type == SignalingMolecule::Drug) {
				if (sm.parents_activating[j]) {
					activating_input_count++;
				}
			} else {
				if (sm.parents_activating[j]) {
					activating_input_count++;
				} else {
					inhibiting_input_count++;
				}
			}
		}

		if (sm.input_mixing_param == SignalingMolecule::Mixing_Undefined && activation_limit == ActivationLimit_LogisticOr) {
			if (activating_input_count > 1) {
				if (!GetVariableIx(varset.get(), std::string("input_mixing_") + sm.name, sm.input_mixing_param)) {
					LOGERROR("\"logistic_or\" activation type is specified, and signaling molecule \"%s\" has more than 1 input but cannot find corresponding mixing parameter", sm.name.c_str());
					return false;
				} else if (sm.input_mixing_param > 50000) {
					LOGERROR("Too many parameters");
					return false;
				}
			} else {
				sm.input_mixing_param = SignalingMolecule::Mixing_Sum;
			}
		}
	}

	// Build a graph structure of the network
	ConstructGraph();

	// Find strongly connected components with Tarjan's algorithm
	component_assignment.resize(signaling_molecules.size());
	num_strongly_connected_components = boost::strong_components(graph, boost::make_iterator_property_map(component_assignment.begin(), boost::get(boost::vertex_index, graph)));
	
	// The strongly connected components are in reverse topological order; reverse the indices to get a topological ordering
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		component_assignment[i] = num_strongly_connected_components - component_assignment[i] - 1;
	}

#if 0
	// Remove mutations - these do not need to be calculated
	std::vector<size_t> remove_ix;
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		const SignalingMolecule& sm = signaling_molecules[i];
		//if (sm.type != SignalingMolecule::Mutation && sm.type != SignalingMolecule::Drug && sm.type != SignalingMolecule::Unknown) {
		if (sm.type == SignalingMolecule::Mutation) {
			size_t component_ix = component_assignment[i];
			remove_ix.push_back(component_ix);
			component_assignment[i] = std::numeric_limits<size_t>::max();
		}
	}
#endif

	if (num_strongly_connected_components >= std::numeric_limits<signaling_index_type>::max() - 1) {
		LOGERROR("Too many strongly connected components!");
		return false;
	}

	// Find the components that contain more than 1 node - these will be solved as systems of equations
	component_size.resize(num_strongly_connected_components);
	component_to_signaling_ix.resize(num_strongly_connected_components);
	std::set<size_t> solve_components;
	for (size_t i = 0; i < num_strongly_connected_components; i++) {
		component_size[i] = std::count(component_assignment.begin(), component_assignment.end(), i);
		if (component_size[i] > 1) {
			solve_components.insert(i);
			component_to_signaling_ix[i] = std::numeric_limits<signaling_index_type>::max();
		} else {
			component_to_signaling_ix[i] = (signaling_index_type)std::distance(component_assignment.begin(), std::find(component_assignment.begin(), component_assignment.end(), i));
		}
	}

	// Allocate structure for system solving
	size_t solve_system_count = 0;
	parallel_data.resize(numthreads);
	for (size_t ti = 0; ti < numthreads; ti++) {
		parallel_data[ti].systems.resize(solve_components.size());
		parallel_data[ti].sobol_sequences.resize(num_strongly_connected_components);

		size_t system_ix = 0;
		for (size_t i = 0; i < num_strongly_connected_components; i++) {
			if (component_size[i] > 1) {
				size_t n = component_size[i];

				NonlinearSystem& system = parallel_data[ti].systems[system_ix];
				system.jac = MatrixReal::Zero(n, n);
				system.vx = VectorReal::Constant(n, 0.5);
				system.fx = VectorReal::Constant(n, 0.0);
				system.deltax = VectorReal::Constant(n, 0.0);
				system.signaling_ix.resize(n);
				size_t j = 0;
				for (size_t k = 0; k < signaling_molecules.size(); k++) {
					if (component_assignment[k] == i) {
						system.signaling_ix[j] = k;
						j++;
					}
				}
				ASSERT(j == n);

				if (n > 4) {
					//parallel_data[ti].systems[system_ix].lu = 
				} else {
					system.inverse = MatrixReal::Zero(n, n);
				}

				parallel_data[ti].sobol_sequences[i] = std::make_unique< boost::random::sobol >(n);

				solve_system_count++;
				system_ix++;
			}
		}

		parallel_data[ti].drug_inhibition.resize(drug_inhibition_counter, std::numeric_limits<Real>::quiet_NaN());
	}

	// Do we  need to solve any systems?
	if (solve_system_count > 0 && activation_limit != ActivationLimit_Logistic) {
		LOGERROR("System contains feedback loop, but the activation limit is not logistic. Feedback loops can only be solved with logistic activation limits.");
		return false;
	}

	// Performance optimizations
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		const SignalingMolecule& sm = signaling_molecules[i];
		if (sm.type == SignalingMolecule::Drug) {
			drug_molecule_indices.push_back(i);
		}
	}

	delete document;
	return true;
}

bool SignalingNetwork::Calculate(size_t threadix, VectorReal& activities, VectorReal& expression, const Real* values)
{
	Precalculate(threadix, activities, expression, values);

	// Calculate network activities
	size_t system_ix = 0;
	for (signaling_index_type i = 0; i < num_strongly_connected_components; i++) {
		signaling_index_type ix = component_to_signaling_ix[i];
		if (ix == std::numeric_limits<signaling_index_type>::max()) {
			// Component contains multiple signaling molecules, we have to solve a system
			NonlinearSystem& system = parallel_data[threadix].systems[system_ix];

			// Always start from the same starting point
			for (size_t i = 0; i < system.signaling_ix.size(); i++) {
				size_t ix = system.signaling_ix[i];
				activities(ix) = 0.5;
			}

			if (!SolveSystem(system, activities, expression, values, threadix)) {
				return false;
			}
			system_ix++;
		} else {
			// Component contains a single signaling molecule, calculate activity directly
			const SignalingMolecule& sm = signaling_molecules[ix];

			if (activities(ix) == activities(ix)) {
				// Activities already set by a model condition or similar (this is only supported outside of feedback loops).
			} else if (sm.type == SignalingMolecule::DrugTransporter) {
				// Nothing to calculate.
			} else {
				Real input, inhibition, maxinput, multinput;
				CalculateActivationInput(ix, activities, values, threadix, input, inhibition, maxinput, multinput);

				if (activation_limit == ActivationLimit_MinMax) {
					Real x = input > 1 ? 1.0 : (input < 0 ? 0 : input);
					activities(ix) = expression_function(input * inhibition, expression(ix), sm.expression_mixing_param, values);
				} else if (activation_limit == ActivationLimit_Logistic) {
					activities(ix) = expression_function(logistic_activation_fixed_bks(input) * inhibition, expression(ix), sm.expression_mixing_param, values);
				} else if (activation_limit == ActivationLimit_LogisticOr) {
					if (sm.input_mixing_param == SignalingMolecule::Mixing_Multiply) {
						activities(ix) = expression_function(logistic_activation_fixed_bks(multinput) * inhibition, expression(ix), sm.expression_mixing_param, values);
					} else if (sm.input_mixing_param == SignalingMolecule::Mixing_Sum) {
						activities(ix) = expression_function(logistic_activation_fixed_bks(input) * inhibition, expression(ix), sm.expression_mixing_param, values);
					} else {
						// Let parameter mix between sum and max
						Real mix = GetValue(values, sm.input_mixing_param);
						activities(ix) = expression_function(logistic_activation_fixed_bks(mix * input + (1.0 - mix) * maxinput) * inhibition, expression(ix), sm.expression_mixing_param, values);
					}
				}
			}
		}
	}

	return true;
}

bool SignalingNetwork::Calculate(size_t threadix, std::vector<VectorReal>& activities, VectorReal& expression, const Real* values)
{
	for (signaling_index_type i = 0; i < num_strongly_connected_components; i++) {
		if (parallel_data[threadix].sobol_sequences[i]) {
			parallel_data[threadix].sobol_sequences[i]->seed();
		}
	}

	// Calculate network activities
	boost::random::uniform_01<Real> unif;
	for (size_t multiroot_i = 0; multiroot_i < activities.size(); multiroot_i++) {
		VectorReal& act = activities[multiroot_i];
		Precalculate(threadix, act, expression, values);

		size_t system_ix = 0;
		for (signaling_index_type i = 0; i < num_strongly_connected_components; i++) {
			signaling_index_type ix = component_to_signaling_ix[i];
			if (ix == std::numeric_limits<signaling_index_type>::max()) {
				// Component contains multiple signaling molecules, we have to solve a system
				NonlinearSystem& system = parallel_data[threadix].systems[system_ix];

				// Get the starting point for this rootfinding iteration
				for (size_t j = 0; j < system.signaling_ix.size(); j++) {
					size_t ix = system.signaling_ix[j];
					act(ix) = unif(*parallel_data[threadix].sobol_sequences[i]);
				}

				if (!SolveSystem(system, act, expression, values, threadix)) {
					return false;
				}
				system_ix++;
			} else {
				// Component contains a single signaling molecule, calculate activity directly
				const SignalingMolecule& sm = signaling_molecules[ix];

				if (act(ix) == act(ix)) {
					// Activities already set by a model condition or similar (this is only supported outside of feedback loops).
				} else if (sm.type == SignalingMolecule::DrugTransporter) {
					// Nothing to calculate.
				} else {
					Real input, inhibition, maxinput, multinput;
					CalculateActivationInput(ix, act, values, threadix, input, inhibition, maxinput, multinput);

					if (activation_limit == ActivationLimit_MinMax) {
						Real x = input > 1 ? 1.0 : (input < 0 ? 0 : input);
						act(ix) = expression_function(input * inhibition, expression(ix), sm.expression_mixing_param, values);
					} else if (activation_limit == ActivationLimit_Logistic) {
						act(ix) = expression_function(logistic_activation_fixed_bks(input) * inhibition, expression(ix), sm.expression_mixing_param, values);
					} else if (activation_limit == ActivationLimit_LogisticOr) {
						if (sm.input_mixing_param == SignalingMolecule::Mixing_Multiply) {
							act(ix) = expression_function(logistic_activation_fixed_bks(multinput) * inhibition, expression(ix), sm.expression_mixing_param, values);
						} else if (sm.input_mixing_param == SignalingMolecule::Mixing_Sum) {
							act(ix) = expression_function(logistic_activation_fixed_bks(input) * inhibition, expression(ix), sm.expression_mixing_param, values);
						} else {
							// Let parameter mix between sum and max
							Real mix = GetValue(values, sm.input_mixing_param);
							act(ix) = expression_function(logistic_activation_fixed_bks(mix * input + (1.0 - mix) * maxinput) * inhibition, expression(ix), sm.expression_mixing_param, values);
						}
					}
				}
			}
		}
	}

	return true;
}

void SignalingNetwork::UpdateActivitiesForDrugSpecies(VectorReal& activities, Real concentration)
{
	for (std::vector<size_t>::const_iterator dmi = drug_molecule_indices.cbegin(); dmi != drug_molecule_indices.cend(); ++dmi) {
		activities(*dmi) = concentration;
	}
}

size_t SignalingNetwork::GetMoleculeCount() const
{
	return signaling_molecules.size();
}

std::string SignalingNetwork::GetMoleculeName(size_t ix) const
{
	if (ix < signaling_molecules.size()) {
		return signaling_molecules[ix].name;
	} else {
		LOGERROR("Invalid molecule index");
		return std::string();
	}
}

size_t SignalingNetwork::GetSignalingMoleculeIxById(const std::string& id)
{
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		if (signaling_molecules[i].id == id) {
			return i;
		}
	}
	return std::numeric_limits<size_t>::max();
}

size_t SignalingNetwork::GetSignalingMoleculeIxByName(const std::string& name)
{
	for (size_t i = 0; i < signaling_molecules.size(); i++) {
		if (signaling_molecules[i].name == name) {
			return i;
		}
	}
	return std::numeric_limits<size_t>::max();
}

void SignalingNetwork::ConstructGraph()
{
	for (signaling_index_type i = 0; i < signaling_molecules.size(); i++) {
		boost::add_vertex(graph);
	}

	for (signaling_index_type i = 0; i < signaling_molecules.size(); i++) {
		const SignalingMolecule& sm = signaling_molecules[i];
		for (signaling_index_type j = 0; j < sm.num_parents; j++) {
			boost::add_edge(sm.parents[j], i, graph);
		}
	}
}

void SignalingNetwork::Precalculate(size_t threadix, VectorReal& activities, VectorReal& expression, const Real* values)
{
	// Handle drug transporter effects
	for (size_t dti = 0; dti < drug_transporters.size(); dti++) {
		const DrugTransporter& dt = drug_transporters[dti];
		for (signaling_index_type k = 0; k < dt.num_targets; k++) {
			if (dt.import[k]) {
				activities(dt.drug_target[k]) = activities(dt.drug_target[k]) * (GetValue(values, dt.strength_varix[k]) * expression(dt.smix) + 1.0 - expression(dt.smix));
			} else {
				activities(dt.drug_target[k]) = activities(dt.drug_target[k]) / (GetValue(values, dt.strength_varix[k]) * expression(dt.smix) + 1.0 - expression(dt.smix));
			}
		}

		// Not actually used in the calculation, but for reference store the expression of the transporter as the activity
		activities(dt.smix) = expression(dt.smix);
	}

	// Precalculate drug inhibition effects
	for (signaling_index_type i = 0; i < signaling_molecules.size(); i++) {
		const SignalingMolecule& sm = signaling_molecules[i];
		for (signaling_index_type j = 0; j < sm.num_parents; j++) {
			signaling_index_type parent_ix = sm.parents[j];
			const SignalingMolecule& parent = signaling_molecules[parent_ix];

			if (parent.type == SignalingMolecule::Drug &&
					(parent.drug_inhibition_type == SignalingMolecule::InhibitActivity ||
					 parent.drug_inhibition_type == SignalingMolecule::InhibitActivityAlterSusceptibility ||
					 parent.drug_inhibition_type == SignalingMolecule::InhibitActivation ||
					 parent.drug_inhibition_type == SignalingMolecule::Activate)) {
				Real signal;
				if (activities[parent_ix] == 0.0) {
					if (sm.parents_activating[j]) {
						signal = 0.0;
					} else {
						signal = 1.0;
					}
				} else {
					Real maxinhib = GetValue(values, sm.param1[j]);
					if (sm.param2[j] == std::numeric_limits<parameter_index_type>::max()) {
						if (sm.parents_activating[j]) {
							signal = activities[parent_ix] * maxinhib;
						} else {
							signal = 1.0 - (activities[parent_ix] * maxinhib);
						}
					} else {
						Real ic50 = GetValue(values, sm.param2[j]);
						Real steepness = bcm3::fastpow10(GetValue(values, sm.param3[j]));
						Real log_concentration = log10(activities[parent_ix]);

						if (sm.parents_activating[j]) {
							// Activating drug
							// maxinhib really means max activation here
							signal = maxinhib - maxinhib / (bcm3::fastpow10(steepness * (log_concentration - ic50)) + 1);
						} else {
							// Inhibiting drug
							signal = (1.0 - maxinhib) + maxinhib / (bcm3::fastpow10(steepness * (log_concentration - ic50)) + 1);
						}
					}
				}
				parallel_data[threadix].drug_inhibition[sm.drug_inhibition_ix[j]] = signal;
			}
		}
	}
}

inline Real SignalingNetwork::CalculateSignalInhibition(size_t i, size_t parent_ix, const VectorReal& activities, const Real* values, size_t threadix) const
{
	const SignalingMolecule& sm = signaling_molecules[i];
	const SignalingMolecule& parent = signaling_molecules[sm.parents[parent_ix]];

	// This is "u_a" in the paper
	Real inhibition = 1.0;

	// Is the parent inhibited by a drug?
	for (signaling_index_type k = 0; k < parent.num_parents; k++) {
		signaling_index_type parent_parent_ix = parent.parents[k];
		const SignalingMolecule& parent_parent = signaling_molecules[parent_parent_ix];
		if (parent_parent.type == SignalingMolecule::Drug &&
			(parent_parent.drug_inhibition_type == SignalingMolecule::InhibitActivity || parent_parent.drug_inhibition_type == SignalingMolecule::InhibitActivityAlterSusceptibility) &&
			activities(parent_parent_ix) > 0 &&
			!parent.parents_activating[k])
		{
			inhibition *= parallel_data[threadix].drug_inhibition[parent.drug_inhibition_ix[k]];
		}
	}

	// Special case -- TODO, this only works correctly if there is only one drug with alter susceptibility
	if (sm.has_susceptibility_altering_parent) {
		for (signaling_index_type k = 0; k < sm.num_parents; k++) {
			const SignalingMolecule& parent2 = signaling_molecules[sm.parents[k]];
			if (parent2.type == SignalingMolecule::Drug &&
				(parent2.drug_inhibition_type == SignalingMolecule::AlterSusceptibility || parent2.drug_inhibition_type == SignalingMolecule::InhibitActivityAlterSusceptibility) &&
				activities(sm.parents[k]) > 0) {
				ASSERT(sm.drug_susceptibility[k] != std::numeric_limits<parameter_index_type>::max());
				inhibition *= GetValue(values, sm.drug_susceptibility[k]);
			}
		}
	}

	return inhibition;
}

inline Real SignalingNetwork::CalculateSignalStrength(size_t i, size_t parent_ix, const VectorReal& activities, const Real* values, size_t threadix) const
{
	const SignalingMolecule& sm = signaling_molecules[i];

	Real signal = GetValue(values, sm.param1[parent_ix]);
	
	if (!sm.parents_activating[parent_ix]) {
		signal = -signal;
	}

	signal *= CalculateSignalInhibition(i, parent_ix, activities, values, threadix);

	return signal;
}

void SignalingNetwork::CalculateActivationInput(size_t i, const VectorReal& activities, const Real* values, size_t threadix, Real& sum, Real& inhibition, Real& maxinput, Real& multinput) const
{
	const SignalingMolecule& sm = signaling_molecules[i];

	// "sum" is the total input, i.e. the result of two sums (one over parent signaling molecules and one over parent mutations) in the paper.
	// "inhibition" is the product over u_b in the paper.

	if (sm.base_parameter_ix == std::numeric_limits<parameter_index_type>::max()) {
		if (sm.num_parents == 0) {
			sum = 1;
		} else {
			sum = 0;
		}
		multinput = 1;
	} else {
		sum = GetValue(values, sm.base_parameter_ix);
		multinput = sum;
	}

	maxinput = sum;
	inhibition = 1.0;
	Real maxinputneg = 0.0;
	Real multinputneg = 0.0;

	// Add signal from all parents which are not inhibiting drugs
	for (signaling_index_type j = 0; j < sm.num_parents; j++) {
		signaling_index_type parent_ix = sm.parents[j];
		const SignalingMolecule& parent = signaling_molecules[parent_ix];

		if (parent.type == SignalingMolecule::Drug) {
			if (sm.parents_activating[j]) {
				Real signal = parallel_data[threadix].drug_inhibition[sm.drug_inhibition_ix[j]];
				sum += signal;
				maxinput = (std::max)(maxinput, signal);
				multinput *= signal;
			} else {
				if (parent.drug_inhibition_type == SignalingMolecule::InhibitActivation || sm.name == "proliferation") {
					Real signal = parallel_data[threadix].drug_inhibition[sm.drug_inhibition_ix[j]];
					inhibition *= signal;
				} else if (parent.drug_inhibition_type == SignalingMolecule::InhibitActivity || parent.drug_inhibition_type == SignalingMolecule::InhibitActivityAlterSusceptibility) {
					// Drugs inhibiting their targets activity only inhibit the downstream signal of their targets, thus they do not affect their immediate target directly
				}
			}
		} else if (parent.type == SignalingMolecule::CompleteLossMutation) {
			if (activities(parent_ix) > 0) {
				sum = 0;
				maxinput = 0;
				multinput = 0;
				break;
			}
		} else {
			Real signal = CalculateSignalStrength(i, j, activities, values, threadix);
			if (sm.param2[j] == std::numeric_limits<parameter_index_type>::max()) {
				// Linear activation
				signal *= activities(parent_ix);
			} else {
				// Non-linear activation
				signal *= logistic_activation(activities(parent_ix), GetValue(values, sm.param3[j]), GetValue(values, sm.param2[j]));
			}
			sum += signal;
			if (sm.parents_activating[j]) {
				maxinput = (std::max)(maxinput, signal);
				multinput *= signal;
			} else {
				maxinputneg = (std::min)(maxinput, signal);
				multinputneg += signal;
			}
		}
	}

	maxinput += maxinputneg;
	multinput += multinputneg;
}

bool SignalingNetwork::SolveSystem(NonlinearSystem& system, VectorReal& activities, VectorReal& expression, const Real* values, size_t threadix) const
{
	static const size_t MAX_NEWTON_ITERATIONS = 20;
	static const Real TOLERANCE = 1e-4;
	
	bool converged = false;

	for (size_t i = 0; i < system.signaling_ix.size(); i++) {
		size_t ix = system.signaling_ix[i];
		system.vx(i) = activities(ix);
	}

	//if (system.signaling_ix.size() == 4) {
	//	LOG("%.8g\t%.8g\t%.8g\t%.8g", system.vx(0), system.vx(1), system.vx(2), system.vx(3));
	//}
	//if (system.signaling_ix.size() == 12) {
	//	LOG("%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g", system.vx(0), system.vx(1), system.vx(2), system.vx(3), system.vx(4), system.vx(5), system.vx(6), system.vx(7), system.vx(8), system.vx(9), system.vx(10), system.vx(11));
	//}
	
	size_t ni = 0;
	while (ni < MAX_NEWTON_ITERATIONS) {
		// Calculate jacobian
		for (size_t i = 0; i < system.signaling_ix.size(); i++) {
			const SignalingMolecule& sm = signaling_molecules[system.signaling_ix[i]];

			Real input, inhibition, maxinput, multinput;
			CalculateActivationInput(system.signaling_ix[i], activities, values, threadix, input, inhibition, maxinput, multinput);

			ASSERT(activation_limit == ActivationLimit_Logistic);

			// F
			system.fx(i) = expression_function(logistic_activation_fixed_bks(input) * inhibition,
				expression(system.signaling_ix[i]), sm.expression_mixing_param, values) - system.vx(i);

			// Jac
			const Real factor = FIXED_K * logistic_deriv_part(input) * inhibition * 
				expression_function(1.0, expression(system.signaling_ix[i]), sm.expression_mixing_param, values);
			for (size_t j = 0; j < system.signaling_ix.size(); j++) {
				if (i == j) {
					system.jac(i,j) = -1.0;
				} else {
					system.jac(i,j) = 0.0;

					for (size_t pixi = 0; pixi < sm.num_parents; pixi++) {
						if (sm.parents[pixi] == system.signaling_ix[j]) {
							Real signal = CalculateSignalStrength(system.signaling_ix[i], pixi, activities, values, threadix);
							system.jac(i,j) = signal * factor;
							break;
						}
					}
				}
			}
		}

		// Solve linear system
		// For size 2-4, use matrix inverse; otherwise LU-decomposition
		switch(system.signaling_ix.size()) {
		case 2:
			Eigen::internal::compute_inverse<MatrixReal, MatrixReal, 2>::run(system.jac, system.inverse);
			system.deltax.noalias() = system.inverse * system.fx;
			break;

		case 3:
			Eigen::internal::compute_inverse<MatrixReal, MatrixReal, 3>::run(system.jac, system.inverse);
			system.deltax.noalias() = system.inverse * system.fx;
			break;

		case 4:
			Eigen::internal::compute_inverse_size4<Eigen::Architecture::SSE, double, MatrixReal, MatrixReal>::run(system.jac, system.inverse);
			system.deltax.noalias() = system.inverse * system.fx;
			break;

		default:
			system.lu.compute_optimized(system.jac);
			system.deltax.noalias() = system.lu.solve(system.fx);
			break;
		}

		// Check whether we've converged and/or diverging
		int conv_count = 0;
		bool diverging = false;
		for (int i = 0; i < system.signaling_ix.size(); i++) {
			Real abschange = fabs(system.deltax(i));
			if (abschange < TOLERANCE) {
				++conv_count;
			} else if (abschange > 0.4) {
				diverging = true;
			}
		}

		// In case of large steps, reduce the change to prevent overshoot
		if (diverging) {
			system.deltax *= 0.5;
		}

		// Calculate new x
		system.vx -= system.deltax;

		for (int i = 0; i < system.signaling_ix.size(); i++) {
			if (system.vx(i) < 0.0) {
				system.vx(i) = 0.0;
			} else if (system.vx(i) > 1.0) {
				system.vx(i) = 1.0;
			}
		}

		//if (system.signaling_ix.size() == 12) {
		//	LOG("%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g", system.vx(0), system.vx(1), system.vx(2), system.vx(3), system.vx(4), system.vx(5), system.vx(6), system.vx(7), system.vx(8), system.vx(9), system.vx(10), system.vx(11));
		//}

		// Update activities vector (this is used in F and Jac calculation)
		for (size_t i = 0; i < system.signaling_ix.size(); i++) {
			size_t ix = system.signaling_ix[i];
			const SignalingMolecule& sm = signaling_molecules[ix];
			activities(ix) = system.vx(i);
		}

		if (conv_count == system.signaling_ix.size()) {
			converged = true;
			break;
		}

		ni++;
	}

	//LOG("Iterations: %u - %s", ni, converged ? "converged" : "NOT converged");
	//if (system.signaling_ix.size() == 12) {
	//	if (converged) {
	//		LOG("%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g", system.vx(0), system.vx(1), system.vx(2), system.vx(3), system.vx(4), system.vx(5), system.vx(6), system.vx(7), system.vx(8), system.vx(9), system.vx(10), system.vx(11));
	//	} else {
	//		LOG("not converged");
	//	}
	//}

	return converged;
}

bool SignalingNetwork::GetVariableIx(const bcm3::VariableSet* varset, const std::string& varstr, parameter_index_type& ix) const
{
	size_t large_ix = varset->GetVariableIndex(varstr);
	if (large_ix == std::numeric_limits<size_t>::max()) {
		return false;
	} else {
		if (large_ix >= std::numeric_limits<parameter_index_type>::max()) {
			LOGERROR("Variable index too large for parameter_index_type");
			return false;
		} else {
			ix = (parameter_index_type)large_ix;
			return true;
		}
	}
}
