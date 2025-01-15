#pragma once

class PharmacokineticModel
{
public:
	PharmacokineticModel();
	~PharmacokineticModel();

	bool SetAbsorption(Real value);
	bool SetExcretion(Real value);
	bool SetElimination(Real value);
	bool SetUsePeripheralCompartment(bool enable);
	bool SetPeripheralForwardRate(Real value);
	bool SetPeripheralBackwardRate(Real value);
	bool SetNumTransitCompartments(size_t value);
	bool SetTransitRate(Real value);

	bool Solve(const VectorReal& treatment_times, const VectorReal& treatment_doses, const VectorReal& observation_timepoints, VectorReal& central_compartment_values, MatrixReal* all_compartments);

private:
	bool UpdateVariable(Real& target, Real value);
	void ConstructMatrix();

	Real absorption;
	Real excretion;
	Real elimination;

	bool use_peripheral_compartment;
	Real peripheral_forward_rate;
	Real peripheral_backward_rate;

	bool use_transit_compartments;
	size_t num_transit_compartments;
	Real transit_rate;

	MatrixReal A;
	MatrixReal tmp1;
	MatrixReal tmp2;
	VectorReal current_y;
	VectorReal tmp_y;

	VectorReal solution_timepoints;
	MatrixReal solution_concentrations;
};
