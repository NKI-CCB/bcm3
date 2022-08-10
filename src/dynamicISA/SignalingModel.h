#pragma once

#include "CVODESolver.h"
#include "SignalingLink.h"
#include "SignalingMolecule.h"

class SignalingModel
{
public:
	SignalingModel();
	~SignalingModel();

	bool Load(const std::string& sbml_fn);
	bool Calculate(const VectorReal& base_param_values,
				   const VectorReal& decay_param_values,
				   const VectorReal& strength_param_values,
				   const VectorReal& inhib_param_values,
				   const VectorReal& expression_mixing_param_values,
				   const VectorReal& fixed_species_values,
				   const VectorReal& expression_levels,
				   Real cell_cycle_duration,
				   VectorReal& activities);
	
	bool Calculate(const VectorReal& base_param_values,
				   const VectorReal& decay_param_values,
				   const VectorReal& strength_param_values,
				   const VectorReal& inhib_param_values,
				   const VectorReal& expression_mixing_param_values,
				   const VectorReal& fixed_species_values,
				   const VectorReal& expression_levels,
				   VectorReal& timepoints,
				   std::vector< std::tuple<size_t, Real, Real> >& treatments,
				   std::vector<VectorReal>& activities);

	inline size_t GetNumMolecules() const { return signaling_molecules.size(); }
	inline size_t GetNumReactions() const { return signaling_links.size(); }
	inline const std::string& GetMoleculeName(size_t i) const { return signaling_molecules[i]->GetName(); }
	std::string GetMoleculeNameById(const std::string& id) const;
	inline const std::string GetStrengthName(size_t i) const { return signaling_links[i]->GetStrengthName(); }
	size_t GetMoleculeIxByName(const std::string& name, bool log_error = true) const;

private:
	bool SetupCalculations(const VectorReal& base_param_values, const VectorReal& fixed_species_values);
	bool CalculateInhibitions();
	bool CalculateDerivative(OdeReal t, const OdeReal* y, OdeReal* dydt, void* user);
	Real TreatmentCallback(OdeReal t, void* user);

	typedef unsigned short index_type;

	struct SolverParams {
		SolverParams();
		const VectorReal* base_param_values;
		const VectorReal* decay_param_values;
		const VectorReal* strength_param_values;
		const VectorReal* inhib_param_values;
		const VectorReal* expression_mixing_param_values;
		const VectorReal* fixed_species_values;
		const VectorReal* expression_levels;
		VectorReal constant_species_values;
		VectorReal inhibitions;

		std::vector< std::tuple<size_t, Real, Real> >* treatments;
		size_t current_discontinuity_ix;
	};

	std::vector< std::unique_ptr<SignalingMolecule> > signaling_molecules;
	std::vector< std::unique_ptr<SignalingLink> > signaling_links;
	CVODESolver solver;
	SolverParams solver_params;

	std::vector< index_type > simulated_species;
};
