#pragma once

namespace bcm3 {

	class Prior;
	class RNG;
	class SampleHistory;

	class Proposal
	{
	public:
		Proposal();
		virtual ~Proposal();

		void Initialize(const std::unique_ptr<SampleHistory>& sample_history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices);
		void GetNewSample(const VectorReal& current_position, VectorReal& new_position, Real& log_mh_ratio, RNG& rng);

		virtual void Update(RNG& rng);
		virtual void NotifyAccepted(bool accepted);

	protected:
		VectorReal Transform(const VectorReal& sample);
		MatrixReal Transform(const MatrixReal& samples);
		VectorReal ReverseTransform(const VectorReal& sample);
		Real ReflectOnBounds(Real x, Real lb, Real ub);
		Real GetPriorVariance(std::shared_ptr<Prior> prior, ptrdiff_t variable_index);

		virtual void InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices) {}
		virtual void GetNewSampleImpl(const VectorReal& current_position, VectorReal& new_position, Real& log_mh_ratio, RNG& rng) = 0;

		enum class ETransforms {
			None,
			Log,
			NegativeLog,
			Logit
		};
		struct Bound {
			Bound();
			Real lower;
			Real upper;
		};

		// Settings
		size_t num_variables;
		bool transform_to_unbounded;
		Real scaling_ema_period;
		Real scaling_learning_rate;
		Real target_acceptance_rate;

		// Runtime variables
		std::vector<ETransforms> variable_transforms;
		std::vector<Bound> variable_bounds;

		Real adaptive_scale;
		Real current_acceptance_rate_ema;
	};

}
