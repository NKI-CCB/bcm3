#pragma once

#include "VariableSet.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/random/sobol.hpp>
#include "EigenPartialPivLUSomewhatSparse.h"

class SignalingNetwork
{
public:
	enum ActivationLimitType
	{
		ActivationLimit_MinMax,
		ActivationLimit_Logistic,
		ActivationLimit_LogisticOr,
		ActivationLimit_Invalid,
	};

	SignalingNetwork(size_t numthreads);
	~SignalingNetwork();

	bool Initialize(const std::string& sbml_file, std::shared_ptr<const bcm3::VariableSet> varset, ActivationLimitType activation_limit, size_t multiroot_solves);
	bool Calculate(size_t threadix, VectorReal& activities, VectorReal& expression, const Real* values);
	bool Calculate(size_t threadix, std::vector<VectorReal>& activities, VectorReal& expression, const Real* values);
	void UpdateActivitiesForDrugSpecies(VectorReal& activities, Real concentration);

	size_t GetMoleculeCount() const;
	std::string GetMoleculeName(size_t ix) const;
	size_t GetSignalingMoleculeIxById(const std::string& id);
	size_t GetSignalingMoleculeIxByName(const std::string& name);

	Real max_expression_function(Real expression, size_t signaling_molecule_ix, const Real* values) const;
	inline size_t GetMultirootSolveCount() const { return multiroot_solves; }
	inline size_t GetNumEvaluationThreads() const { return numthreads; }

private:
	typedef unsigned char signaling_index_type;
	typedef unsigned short parameter_index_type;
	static const int MAX_PARENTS = 16;

	struct SignalingMolecule {
		SignalingMolecule();

		enum EType {
			Protein,
			mRNA,
			SmallMolecule,
			Mutation,
			CompleteLossMutation,
			Drug,
			Phenotype,
			Unknown,
			DrugTransporter,
			InvalidType,
		};
		enum EDrugInhibitionType {
			InhibitActivity,
			InhibitActivityAlterSusceptibility,
			InhibitActivation,
			AlterSusceptibility,
			Activate,
			NotADrug,
			InvalidDrugInhibitionType,
		};
		enum EInputMixingType {
			Mixing_Multiply = 65533,
			Mixing_Sum = 65534,
			Mixing_Undefined = 65535,
		};

		std::string id;
		std::string name;
		
		char type;
		char drug_inhibition_type;
		parameter_index_type base_parameter_ix;

		signaling_index_type num_parents;
		signaling_index_type parents[MAX_PARENTS];
		parameter_index_type param1[MAX_PARENTS]; // linear: strength	-- nonlinear: strength		-- drug: maxinhib
		parameter_index_type param2[MAX_PARENTS]; // linear: nothing	-- nonlinear: inflection	-- drug: ic50
		parameter_index_type param3[MAX_PARENTS]; // linear: nothing	-- nonlinear: steepness		-- drug: steepness
		signaling_index_type drug_inhibition_ix[MAX_PARENTS];
		unsigned char parents_activating[MAX_PARENTS];
		parameter_index_type drug_susceptibility[MAX_PARENTS];
		parameter_index_type input_mixing_param;
		parameter_index_type expression_mixing_param;
		bool has_susceptibility_altering_parent; // performance optimization
	};

	struct DrugTransporter {
		DrugTransporter();

		std::string id;
		std::string name;

		signaling_index_type smix;
		signaling_index_type num_targets;
		signaling_index_type drug_target[MAX_PARENTS];
		parameter_index_type strength_varix[MAX_PARENTS];
		unsigned char import[MAX_PARENTS];
	};

	struct NonlinearSystem {
		PartialPivLUExtended lu;
		VectorReal vx, fx, deltax;
		MatrixReal jac;
		MatrixReal inverse;
		std::vector<size_t> signaling_ix;
	};
	
	struct ParallelData {
		std::vector<NonlinearSystem> systems;
		std::vector<Real> drug_inhibition;
		std::vector< std::unique_ptr<boost::random::sobol> > sobol_sequences;
	};

	void ConstructGraph();
	void Precalculate(size_t threadix, VectorReal& activities, VectorReal& expression, const Real* values);
	void CalculateActivationInput(size_t i, const VectorReal& activities, const Real* values, size_t threadix, Real& sum, Real& inhibition, Real& maxinput, Real& multinput) const;
	inline Real CalculateSignalInhibition(size_t i, size_t parent_ix, const VectorReal& activities, const Real* values, size_t threadix) const;
	inline Real CalculateSignalStrength(size_t i, size_t parent_ix, const VectorReal& activities, const Real* values, size_t threadix) const;
	bool SolveSystem(NonlinearSystem& system, VectorReal& activities, VectorReal& expression, const Real* values, size_t threadix) const;
	bool GetVariableIx(const bcm3::VariableSet* varset, const std::string& varstr, parameter_index_type& ix) const;
	inline Real expression_function(Real activity, Real expression, SignalingNetwork::parameter_index_type expression_mixing_param, const Real* values) const;
	
	inline Real GetValue(const Real* values, parameter_index_type varix) const {
		ASSERT(varix < VarSet->GetNumVariables());
		return values[varix];
	}

	size_t numthreads;
	std::vector<ParallelData> parallel_data;
	
	std::shared_ptr<const bcm3::VariableSet> VarSet;
	ActivationLimitType activation_limit;
	size_t multiroot_solves;

	std::vector<SignalingMolecule> signaling_molecules;
	std::vector<DrugTransporter> drug_transporters;

	boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> graph;
	size_t num_strongly_connected_components;
	std::vector<size_t> component_assignment;
	std::vector<size_t> component_size;
	std::vector<signaling_index_type> component_to_signaling_ix;

	// Optimizations
	std::vector<size_t> drug_molecule_indices;
};
