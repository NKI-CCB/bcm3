#pragma once

namespace bcm3 {

	class Prior;
	class RNG;
	class SampleHistory;
	class SampleHistoryClustering;

	class Proposal
	{
	public:
		Proposal();
		virtual ~Proposal();

		bool Initialize(const SampleHistory& sample_history, const std::shared_ptr<const SampleHistoryClustering> sample_history_clustering,
			size_t max_history_samples, bool transform_to_unbounded, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng,
			const std::string& tmpfilename, bool log_info);

		virtual void GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng) = 0;
		virtual Real CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment) = 0;

		virtual bool UsesClustering();
		virtual void Update(RNG& rng);
		virtual void NotifyAccepted(bool accepted);
		virtual void LogInfo() const;
		virtual void WriteToFile(const std::string& fn, std::string adaptation_group, std::vector<ptrdiff_t>& variable_indices);

	protected:
#if 0
		VectorReal Transform(const VectorReal& sample);
		MatrixReal Transform(const MatrixReal& samples);
		VectorReal ReverseTransform(const VectorReal& sample);
#endif
		Real ReflectOnBounds(Real x, Real lb, Real ub);
		Real GetPriorVariance(std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, ptrdiff_t i);

		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);

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
		std::string tmpfilename;

		// Runtime variables
		std::shared_ptr<const SampleHistoryClustering> sample_history_clustering;
		std::vector<ETransforms> variable_transforms;
		std::vector<Bound> variable_bounds;

		Real adaptive_scale;
		Real current_acceptance_rate_ema;
	};

}
