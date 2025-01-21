#include "PharmacokineticModel.h"

#include "unsupported/Eigen/MatrixFunctions"

PharmacokineticModel::PharmacokineticModel()
	: bioavailability(1.0)
	, absorption(std::numeric_limits<Real>::quiet_NaN())
	, excretion(std::numeric_limits<Real>::quiet_NaN())
	, elimination(std::numeric_limits<Real>::quiet_NaN())
	, use_peripheral_compartment(false)
	, peripheral_forward_rate(std::numeric_limits<Real>::quiet_NaN())
	, peripheral_backward_rate(std::numeric_limits<Real>::quiet_NaN())
	, use_transit_compartments(false)
	, num_transit_compartments(0)
	, transit_rate(std::numeric_limits<Real>::quiet_NaN())
	, use_biphasic_abosprtion(false)
	, direct_absorption_rate(std::numeric_limits<Real>::quiet_NaN())
	, fraction_direct(std::numeric_limits<Real>::quiet_NaN())
	, use_metabolite(false)
	, metabolite_conversion_rate(std::numeric_limits<Real>::quiet_NaN())
	, metabolite_elimination(std::numeric_limits<Real>::quiet_NaN())
{
}

PharmacokineticModel::~PharmacokineticModel()
{
}

bool PharmacokineticModel::SetBioavailability(Real value)
{
	return UpdateVariable(this->bioavailability, value);
}

bool PharmacokineticModel::SetAbsorption(Real value)
{
	return UpdateVariable(this->absorption, value);
}

bool PharmacokineticModel::SetExcretion(Real value)
{
	return UpdateVariable(this->excretion, value);
}

bool PharmacokineticModel::SetElimination(Real value)
{
	return UpdateVariable(this->elimination, value);
}

bool PharmacokineticModel::SetUsePeripheralCompartment(bool enable)
{
	use_peripheral_compartment = enable;
	return true;
}

bool PharmacokineticModel::SetPeripheralForwardRate(Real value)
{
	return UpdateVariable(this->peripheral_forward_rate, value);
}

bool PharmacokineticModel::SetPeripheralBackwardRate(Real value)
{
	return UpdateVariable(this->peripheral_backward_rate, value);
}

bool PharmacokineticModel::SetNumTransitCompartments(size_t value)
{
	if (this->num_transit_compartments != value) {
		if (value > 0) {
			use_transit_compartments = true;
		} else {
			use_transit_compartments = false;
		}
		this->num_transit_compartments = value;
		return true;
	} else {
		return false;
	}
}

bool PharmacokineticModel::SetTransitRate(Real value)
{
	return UpdateVariable(this->transit_rate, value);
}

bool PharmacokineticModel::SetUseBiphasicAbsorption(bool enable)
{
	use_biphasic_abosprtion = enable;
	return true;
}

bool PharmacokineticModel::SetDirectAbsorptionRate(Real value)
{
	return UpdateVariable(this->direct_absorption_rate, value);
}

bool PharmacokineticModel::SetFractionDirect(Real value)
{
	return UpdateVariable(this->fraction_direct, value);
}

bool PharmacokineticModel::SetUseMetabolite(bool enable)
{
	use_metabolite = enable;
	return true;
}

bool PharmacokineticModel::SetMetaboliteConversionRate(Real value)
{
	return UpdateVariable(this->metabolite_conversion_rate, value);
}

bool PharmacokineticModel::SetMetaboliteElimination(Real value)
{
	return UpdateVariable(this->metabolite_elimination, value);
}

bool PharmacokineticModel::Solve(const VectorReal& treatment_times, const VectorReal& treatment_doses, const VectorReal& observation_timepoints, VectorReal& central_compartment_values, MatrixReal* all_compartments)
{
	ASSERT(treatment_doses.size() == treatment_times.size());
	ASSERT(central_compartment_values.size() == observation_timepoints.size());

	ConstructMatrix();

	solution_timepoints = treatment_times;
	solution_concentrations.resize(treatment_times.size(), A.rows());

	if (all_compartments) {
		all_compartments->setConstant(observation_timepoints.size(), A.rows(), std::numeric_limits<Real>::quiet_NaN());
	}

	current_y.setZero(A.rows());

	Real simulate_until = observation_timepoints.tail(1)(0);

	ptrdiff_t tti = 0;
	ptrdiff_t oti = 0;
	Real current_t = 0;
	while (tti < treatment_times.size() && current_t < simulate_until) {
		// Calculate target time for this timestep
		Real target_t;
		if (tti < treatment_times.size() - 1) {
			target_t = treatment_times[tti + 1];
		} else {
			target_t = simulate_until;
		}
		current_y(0) += treatment_doses(tti) * bioavailability;

		// If there are any observation timepoints between now and the next treatment timepoint (or end of simulation), evaluate them
		while (oti < observation_timepoints.size() && observation_timepoints(oti) <= target_t) {
			Real offset_t = observation_timepoints(oti) - current_t;
			tmp1.noalias() = A * offset_t;
			tmp2.noalias() = tmp1.exp();

			tmp_y = tmp2 * current_y;
			central_compartment_values(oti) = tmp_y(1);

			if (all_compartments) {
				all_compartments->row(oti) = tmp_y;
			}

			oti++;
		}

		// Advance the trajectory to the next treatment timepoint
		Real dt = target_t - current_t;
		tmp1.noalias() = A * dt;
		tmp2.noalias()  = tmp1.exp();

		solution_concentrations.row(tti) = tmp2 * current_y;

		for (int i = 0; i < A.rows(); i++) {
			if (std::isnan(solution_concentrations(tti, i))) {
				return false;
			}
		}

		current_t = target_t;
		current_y = solution_concentrations.row(tti);
		tti++;
	}

	return true;
}

bool PharmacokineticModel::UpdateVariable(Real& target, Real value)
{
	if (target != value) {
		target = value;
		return true;
	} else {
		return false;
	}
}

void PharmacokineticModel::ConstructMatrix()
{
	size_t num_compartments = 2;
	size_t metabolite_compartment_ix = std::numeric_limits<size_t>::max();
	size_t first_transit_compartment_ix = 0;
	if (use_peripheral_compartment) {
		num_compartments++;
	}
	if (use_metabolite) {
		metabolite_compartment_ix = num_compartments;
		num_compartments++;
	}
	if (use_transit_compartments) {
		first_transit_compartment_ix = num_compartments;
		num_compartments += num_transit_compartments;
	}

	A.setZero(num_compartments, num_compartments);

	A(0, 0) -= excretion;

	Real absorption_indirect;
	if (use_biphasic_abosprtion) {
		absorption_indirect = (1 - fraction_direct) * absorption;
	} else {
		absorption_indirect = absorption;
	}
	A(0, 0) -= absorption_indirect;
	
	if (use_transit_compartments) {
		A(first_transit_compartment_ix, 0) += absorption_indirect;
		if (num_transit_compartments > 2) {
			for (size_t i = 0; i < num_transit_compartments - 1; i++) {
				A(first_transit_compartment_ix + i,     first_transit_compartment_ix + i) -= transit_rate;
				A(first_transit_compartment_ix + i + 1, first_transit_compartment_ix + i) += transit_rate;
			}
		}
		A(first_transit_compartment_ix + num_transit_compartments - 1, first_transit_compartment_ix + num_transit_compartments - 1) = -transit_rate;
		A(1, first_transit_compartment_ix + num_transit_compartments - 1) += transit_rate;
	} else {
		A(1, 0) += absorption_indirect;
	}

	if (use_peripheral_compartment) {
		A(1, 1) -= peripheral_forward_rate;
		A(2, 1) += peripheral_forward_rate;
		A(1, 2) += peripheral_backward_rate;
		A(2, 2) -= peripheral_backward_rate;
	}

	if (use_biphasic_abosprtion) {
		A(0, 0) -= fraction_direct * direct_absorption_rate;
		A(1, 0) += fraction_direct * direct_absorption_rate;
	}

	if (use_metabolite) {
		A(1, 1) -= metabolite_conversion_rate;
		A(metabolite_compartment_ix, 1) += metabolite_conversion_rate;
		A(metabolite_compartment_ix, metabolite_compartment_ix) -= metabolite_elimination;
	}

	A(1, 1) -= elimination;

	tmp1.resizeLike(A);
	tmp2.resizeLike(A);
}
