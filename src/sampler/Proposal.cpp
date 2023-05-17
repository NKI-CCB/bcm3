#include "Utils.h"
#include "Prior.h"
#include "Proposal.h"
#include "RNG.h"
#include "SampleHistory.h"
#include "SampleHistoryClustering.h"

namespace bcm3 {

	// The variable bounds are inclusive so values can be exactly at the bounds
	// This could give infinite transformed values, so we extend the logit scale slightly beyond the bounds in the forward transform,
	// and reflect on the bounds in the reverse transform
	static const Real logit_transform_margin = 1e-5;

	Proposal::Bound::Bound()
		: lower(std::numeric_limits<Real>::quiet_NaN())
		, upper(std::numeric_limits<Real>::quiet_NaN())
	{
	}

	Proposal::Proposal()
		: num_variables(0)
		, transform_to_unbounded(false)
		, scaling_ema_period(1000)
		, scaling_learning_rate(0.05)
		, target_acceptance_rate(0.234)
		, t_dof(0.0)
		, adaptive_scale(1.0)
		, current_acceptance_rate_ema(0.23)
	{
	}

	Proposal::~Proposal()
	{
	}

	bool Proposal::Initialize(const SampleHistory& sample_history, const std::shared_ptr<SampleHistoryClustering> sample_history_clustering,
		size_t max_history_samples, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, Real t_dof,
		const std::string& tmpfilename, bool log_info)
	{
		num_variables = variable_indices.size();
		this->sample_history_clustering = sample_history_clustering;
		this->tmpfilename = tmpfilename;
		this->t_dof = t_dof;

		if (num_variables == 1) {
			target_acceptance_rate = 0.44;
		} else if (num_variables == 2) {
			target_acceptance_rate = 0.35;
		} else if (num_variables == 3) {
			target_acceptance_rate = 0.3;
		} else {
			target_acceptance_rate = 0.234;
		}

		variable_bounds.resize(num_variables);
		this->transform_to_unbounded = transform_to_unbounded;
		if (transform_to_unbounded) {
			variable_transforms.resize(num_variables);
		}

		for (size_t i = 0; i < variable_indices.size(); i++) {
			variable_bounds[i].lower = prior->GetLowerBound(variable_indices[i]);
			variable_bounds[i].upper = prior->GetUpperBound(variable_indices[i]);

			if (transform_to_unbounded) {
				if (variable_bounds[i].lower == -std::numeric_limits<Real>::infinity() && variable_bounds[i].upper == std::numeric_limits<Real>::infinity()) {
					variable_transforms[i] = ETransforms::None;
				} else if (variable_bounds[i].lower == -std::numeric_limits<Real>::infinity() && variable_bounds[i].upper != std::numeric_limits<Real>::infinity()) {
					variable_transforms[i] = ETransforms::NegativeLog;
				} else if (variable_bounds[i].lower != -std::numeric_limits<Real>::infinity() && variable_bounds[i].upper == std::numeric_limits<Real>::infinity()) {
					variable_transforms[i] = ETransforms::Log;
				} else {
					variable_transforms[i] = ETransforms::Logit;
				}
			}
		}

		MatrixReal history = sample_history.GetHistory(variable_indices);
		if (log_info) {
			LOG("  Proposal adaptation - history samples: %zd", history.rows());
		}

		// If necessary, reduce the history size to the maximum number of samples specified in the configuration
		if (history.rows() > max_history_samples) {
			if (log_info) {
				LOG("Proposal adaptation - downsampling to %zu samples", max_history_samples);
			}

			std::vector<ptrdiff_t> use_sample_ix;

			// First subsample to get at most n*2-1 samples
			ptrdiff_t subsample = history.rows() / max_history_samples;
			if (subsample > 1) {
				use_sample_ix.resize(history.rows() / subsample);
				for (ptrdiff_t i = 0; i < use_sample_ix.size(); i++) {
					use_sample_ix[i] = i * subsample;
				}
			} else {
				use_sample_ix.resize(history.rows());
				for (ptrdiff_t i = 0; i < use_sample_ix.size(); i++) {
					use_sample_ix[i] = i;
				}
			}

			// Then discard samples at random until we have the required amount
			while (use_sample_ix.size() > max_history_samples) {
				unsigned int drop_sample = rng.GetUnsignedInt(use_sample_ix.size() - 1);
				use_sample_ix.erase(use_sample_ix.begin() + drop_sample);
			}

			MatrixReal selected = history(use_sample_ix, Eigen::all);
			history = selected;

			if (sample_history_clustering) {
				sample_history_clustering->AssignAllHistorySamples(use_sample_ix, history);
			}
		} else {
			if (sample_history_clustering) {
				std::vector<ptrdiff_t> use_sample_ix(history.rows());
				for (ptrdiff_t i = 0; i < use_sample_ix.size(); i++) {
					use_sample_ix[i] = i;
				}
				sample_history_clustering->AssignAllHistorySamples(use_sample_ix, history);
			}
		}

#if 0
		if (transform_to_unbounded) {
			history = Transform(history);
		}
#else
		ASSERT(!transform_to_unbounded);
#endif

		return InitializeImpl(history, prior, variable_indices, rng, log_info);
	}

#if 0
	// Transformation don't seem to add much - deprecate for now
	void Proposal::GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng)
	{
		ASSERT(current_position.size() == num_variables);

		if (transform_to_unbounded) {
			VectorReal transformed = Transform(current_position);
			GetNewSampleImpl(transformed, history_cluster_assignment, new_position, rng);
			new_position = ReverseTransform(new_position);

			for (ptrdiff_t i = 0; i < num_variables; i++) {
				switch (variable_transforms[i]) {
				case ETransforms::None:
					break;

				case ETransforms::Log:
					{
						// Derivative of log(x - l) = 1 / (x-l)
						Real l = variable_bounds[i].lower;
						log_mh_ratio += (-log(current_position[i] - l)) - (-log(new_position[i] - l));
					}
					break;

				case ETransforms::NegativeLog:
					{
						// Derivative of log(u - x) = 1 / (x-l)
						Real u = variable_bounds[i].upper;
						log_mh_ratio += (-log(u - current_position[i])) - (-log(u - new_position[i]));
					}
					break;

				case ETransforms::Logit:
					{
						Real l = variable_bounds[i].lower;
						Real u = variable_bounds[i].upper;
						Real margin = logit_transform_margin * (u - l);
						log_mh_ratio += log(dlogit_scale(current_position[i], l-margin, u+margin)) - log(dlogit_scale(new_position[i], l-margin, u+margin));
					}
					break;

				default:
					LOGERROR("Invalid variable transform");
				}
			}
		} else {
			GetNewSampleImpl(current_position, history_cluster_assignment, new_position, rng);
			for (ptrdiff_t i = 0; i < num_variables; i++) {
				new_position(i) = ReflectOnBounds(new_position(i), variable_bounds[i].lower, variable_bounds[i].upper);
			}
		}
	}
#endif

	bool Proposal::UsesClustering()
	{
		return false;
	}

	void Proposal::Update(RNG& rng)
	{
		Real learn_rate = 1.0 + rng.GetReal() * scaling_learning_rate;
		if (current_acceptance_rate_ema < 0.952381 * target_acceptance_rate) {
			adaptive_scale /= learn_rate;
			adaptive_scale = std::max(adaptive_scale, (Real)1e-4);
		} else if (current_acceptance_rate_ema > 1.05 * target_acceptance_rate) {
			adaptive_scale *= learn_rate;
			adaptive_scale = std::min(adaptive_scale, (Real)10.0);
		}
	}

	void Proposal::NotifyAccepted(bool accepted)
	{
		const Real ema_alpha = 2.0 / (scaling_ema_period + 1);

		if (accepted) {
			current_acceptance_rate_ema += (1.0 - current_acceptance_rate_ema) * ema_alpha;
		} else {
			current_acceptance_rate_ema += (0.0 - current_acceptance_rate_ema) * ema_alpha;
		}
	}

	void Proposal::LogInfo() const
	{
	}

	void Proposal::WriteToFile(const std::string& fn, std::string adaptation_group, std::vector<ptrdiff_t>& variable_indices)
	{
	}

#if 0
	VectorReal Proposal::Transform(const VectorReal& sample)
	{
		ASSERT(sample.size() == variable_transforms.size());

		if (transform_to_unbounded) {
			VectorReal transformed(sample.size());
			for (ptrdiff_t i = 0; i < sample.size(); i++) {
				switch (variable_transforms[i]) {
				case ETransforms::None:
					transformed[i] = sample[i];
					break;

				case ETransforms::Log:
					transformed[i] = log(sample[i] - variable_bounds[i].lower);
					break;

				case ETransforms::NegativeLog:
					transformed[i] = log(variable_bounds[i].upper - sample[i]);
					break;

				case ETransforms::Logit:
					{
						Real range = variable_bounds[i].upper - variable_bounds[i].lower;
						Real margin = logit_transform_margin * range;
						Real p = (sample[i] - (variable_bounds[i].lower - margin)) / (range + 2 * margin);
						transformed[i] = log(p / (1.0 - p));
					}
					break;

				default:
					LOGERROR("Invalid variable transform");
				}

				ASSERT(!std::isnan(sample(i)));
			}

			return transformed;
		} else {
			return sample;
		}
	}

	MatrixReal Proposal::Transform(const MatrixReal& samples)
	{
		if (samples.rows() == 0) {
			return samples;
		}

		ASSERT(samples.cols() == variable_transforms.size());

		if (transform_to_unbounded) {
			MatrixReal transformed(samples.rows(), samples.cols());
			for (ptrdiff_t i = 0; i < samples.cols(); i++) {
				switch (variable_transforms[i]) {
				case ETransforms::None:
					transformed.col(i) = samples.col(i);
					break;

				case ETransforms::Log:
					transformed.col(i) = (samples.col(i).array() - variable_bounds[i].lower).log();
					break;

				case ETransforms::NegativeLog:
					transformed.col(i) = (variable_bounds[i].upper - samples.col(i).array()).log();
					break;

				case ETransforms::Logit:
				{
					Real range = variable_bounds[i].upper - variable_bounds[i].lower;
					Real margin = logit_transform_margin * range;
					VectorReal p = (samples.col(i).array() - (variable_bounds[i].lower - margin)) / (range + 2 * margin);
					transformed.col(i) = (p.array() / (1.0 - p.array())).log();
				}
				break;

				default:
					LOGERROR("Invalid variable transform");
				}

				for (ptrdiff_t j = 0; j < samples.rows(); j++) {
					ASSERT(!std::isnan(transformed(j, i)));
				}
			}
			return transformed;
		} else {
			return samples;
		}
	}

	VectorReal Proposal::ReverseTransform(const VectorReal& sample)
	{
		ASSERT(sample.size() == variable_transforms.size());

		if (transform_to_unbounded) {
			VectorReal reverse_transformed(sample.size());
			for (ptrdiff_t i = 0; i < sample.size(); i++) {
				switch (variable_transforms[i]) {
				case ETransforms::None:
					reverse_transformed(i) = sample(i);
					break;

				case ETransforms::Log:
					reverse_transformed(i) = variable_bounds[i].lower + exp(sample(i));
					break;

				case ETransforms::NegativeLog:
					reverse_transformed(i) = variable_bounds[i].upper - exp(sample(i));
					break;

				case ETransforms::Logit:
					{
						Real range = variable_bounds[i].upper - variable_bounds[i].lower;
						Real margin = logit_transform_margin * range;
						Real p = 1.0 / (1.0 + exp(-sample(i)));
						Real x = (variable_bounds[i].lower - margin) + p * (range - 2 * margin);
						reverse_transformed(i) = ReflectOnBounds(x, variable_bounds[i].lower, variable_bounds[i].upper);
					}
					break;

				default:
					LOGERROR("Invalid variable transform");
				}

				ASSERT(!std::isnan(reverse_transformed(i)));
			}
			return reverse_transformed;
		} else {
			return sample;
		}
	}
#endif

	Real Proposal::ReflectOnBounds(Real x, Real lb, Real ub)
	{
		while (1) {
			if (x < lb) {
				x = lb + (lb - x);
			} else if (x > ub) {
				x = ub - (x - ub);
			} else {
				break;
			}
		}
		return x;
	}

	Real Proposal::GetPriorVariance(std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, ptrdiff_t i)
	{
		Real var = std::numeric_limits<Real>::quiet_NaN();
		if (transform_to_unbounded) {
			switch (variable_transforms[i]) {
			case ETransforms::None:
				prior->EvaluateMarginalVariance(variable_indices[i], var);
				break;

			case ETransforms::Log:
				ASSERT(false);
				var = 1.0;
				break;

			case ETransforms::NegativeLog:
				ASSERT(false);
				var = 1.0;
				break;

			case ETransforms::Logit:
				// Assume it's a uniform prior distribution for now; in which case the variance of the logit transformed variable is the variance of the logistic distribution
				var = M_PI * M_PI / 3.0;
				break;

			default:
				LOGERROR("Invalid variable transform");
			}
		} else {
			prior->EvaluateMarginalVariance(variable_indices[i], var);
		}
		return var;
	}

	bool Proposal::InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info)
	{
		return true;
	}
}
