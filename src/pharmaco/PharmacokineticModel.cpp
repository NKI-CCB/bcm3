#include "PharmacokineticModel.h"

#include "unsupported/Eigen/MatrixFunctions"

PharmacokineticModel::PharmacokineticModel()
	: absorption(std::numeric_limits<Real>::quiet_NaN())
	, excretion(std::numeric_limits<Real>::quiet_NaN())
	, elimination(std::numeric_limits<Real>::quiet_NaN())
{

}

PharmacokineticModel::~PharmacokineticModel()
{

}

bool PharmacokineticModel::SetAbsorption(Real value)
{
	if (this->absorption != value) {
		this->absorption = value;
		return true;
	} else {
		return false;
	}
}

bool PharmacokineticModel::SetExcretion(Real value)
{
	if (this->excretion != value) {
		this->excretion = value;
		return true;
	} else {
		return false;
	}
}

bool PharmacokineticModel::SetElimination(Real value)
{
	if (this->elimination != value) {
		this->elimination = value;
		return true;
	} else {
		return false;
	}
}

bool PharmacokineticModel::Solve(const VectorReal& treatment_times, const VectorReal& treatment_doses, const VectorReal& observation_timepoints, VectorReal& central_compartment_values)
{
	ASSERT(treatment_doses.size() == treatment_times.size());
	ASSERT(central_compartment_values.size() == observation_timepoints.size());

	ConstructMatrix();

	solution_timepoints = treatment_times;
	solution_concentrations.resize(treatment_times.size(), 2);

	current_y.setZero(2);

	ptrdiff_t oti = 0;
	Real current_t = 0;
	for (ptrdiff_t i = 0; i < treatment_times.size(); i++) {
		// If there are any observation timepoints between now and the next treatment timepoint (or end of simulation), evalute them
		Real evaluate_observation_until;
		if (i == treatment_times.size() - 1) {
			evaluate_observation_until = std::numeric_limits<Real>::max();
		} else {
			evaluate_observation_until = treatment_times(i + 1);
		}
		while (oti < observation_timepoints.size() && observation_timepoints(oti) < evaluate_observation_until) {
			Real offset_t = observation_timepoints(oti) - current_t;
			tmp1.noalias() = A * offset_t;
			tmp2.noalias() = tmp1.exp();

			tmp_y = tmp2 * current_y;
			central_compartment_values(oti) = tmp_y(1);
			oti++;
		}

		// Advance the trajectory to the next treatment timepoint
		Real dt;
		if (i == 0) {
			dt = treatment_times[0];
		} else {
			dt = treatment_times[i] - treatment_times[i - 1];
		}

		tmp1.noalias() = A * dt;
		tmp2.noalias()  = tmp1.exp();

		current_y(0) += treatment_doses(i);

		solution_concentrations.row(i) = tmp2 * current_y;

		current_t = treatment_times(i);
		current_y = solution_concentrations.row(i);
	}

	return true;
}

void PharmacokineticModel::Evaluate(Real t, Real& central_y)
{

}

void PharmacokineticModel::Evaluate(Real t, VectorReal& y)
{

}

void PharmacokineticModel::ConstructMatrix()
{
	A.resize(2, 2);
	A(0, 0) = -(absorption + excretion);
	A(0, 1) = 0;
	A(1, 0) = absorption;
	A(1, 1) = -elimination;

	tmp1.resizeLike(A);
	tmp2.resizeLike(A);
}
