#include "Utils.h"
#include "Prior.h"
#include "ProposalGlobalCovariance.h"
#include "RNG.h"
#include "SummaryStats.h"

namespace bcm3 {

	ProposalGlobalCovariance::ProposalGlobalCovariance()
		: logC(std::numeric_limits<Real>::quiet_NaN())
		, t_dof(0.0)
	{
	}

	ProposalGlobalCovariance::~ProposalGlobalCovariance()
	{
	}

	void ProposalGlobalCovariance::GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng)
	{
		ASSERT(current_position.size() == num_variables);

		Real t_scale;
		if (t_dof > 0.0) {
			Real w = rng.GetGamma(0.5 * t_dof, 0.5 * t_dof);
			t_scale = bcm3::rsqrt(w);
		} else {
			t_scale = 1.0;
		}

		VectorReal x = rng.GetMultivariateUnitNormal(num_variables);
		x = covariance_decomp * x;
		x *= t_scale * adaptive_scale;
		new_position = current_position + x;

		for (ptrdiff_t i = 0; i < num_variables; i++) {
			new_position(i) = ReflectOnBounds(new_position(i), variable_bounds[i].lower, variable_bounds[i].upper);
		}
	}

	Real ProposalGlobalCovariance::CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment)
	{
		// Proposal is symmetric
		return 0.0;
	}

	void ProposalGlobalCovariance::LogInfo() const
	{
		LOG("  Global covariance; scale=%8.5f, condition number=%6g", adaptive_scale, 1.0 / covariance_llt.rcond());
	}

	bool ProposalGlobalCovariance::InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info)
	{
		if (history.rows() < 2) {
			covariance = MatrixReal::Identity(num_variables, num_variables);
			for (size_t j = 0; j < variable_indices.size(); j++) {
				covariance(j,j) = GetPriorVariance(prior, variable_indices, j);
			}
		} else {
			covariance = cov(history);

#if 0
			std::stringstream str;
			str << covariance;
			LOG("Covariance:");
			LOG("\n%s", str.str().c_str());
#endif

			for (ptrdiff_t i = 0; i < variable_indices.size(); i++) {
				ASSERT(!std::isnan(covariance(i,i)));
			}

			// Make sure the diagonal is always at least a small positive value
			for (size_t j = 0; j < variable_indices.size(); j++) {
				Real prior_var = GetPriorVariance(prior, variable_indices, j);
				covariance(j, j) = std::max(covariance(j, j), 1e-6 * prior_var);
			}
		}

		// Cholesky decomposition
		covariance_llt = covariance.llt();
		covariance_decomp = covariance_llt.matrixL();

		// Precalculate the normalization constant
		Real det = 0.0;
		for (size_t j = 0; j < variable_indices.size(); j++) {
			det += log(covariance_decomp(j, j));
		}
		logC = -det - 0.5 * variable_indices.size() * log(2.0 * M_PI);

		return true;
	}

}
