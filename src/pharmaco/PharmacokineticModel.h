#pragma once

class PharmacokineticModel
{
public:
	PharmacokineticModel();
	~PharmacokineticModel();

	bool SetAbsorption(Real value);
	bool SetExcretion(Real value);
	bool SetElimination(Real value);
	bool SetIntercompartmentalForward(Real value);
	bool SetIntercompartmentalBackward(Real value);
	bool SetNumTransitCompartments(size_t value);
	bool SetMeanTransitTime(Real value);

	bool Solve(const VectorReal& treatment_times, const VectorReal& treatment_doses, const VectorReal& observation_timepoints, VectorReal& central_compartment_values);

	void Evaluate(Real t, Real& central_y);
	void Evaluate(Real t, VectorReal& y);

private:
	void ConstructMatrix();

	Real absorption;
	Real excretion;
	Real elimination;

	MatrixReal A;
	MatrixReal tmp1;
	MatrixReal tmp2;
	VectorReal current_y;
	VectorReal tmp_y;

	VectorReal solution_timepoints;
	MatrixReal solution_concentrations;
};
