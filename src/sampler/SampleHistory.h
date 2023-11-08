#pragma once

namespace bcm3 {

	class SampleHistory
	{
	public:
		SampleHistory();
		~SampleHistory();

		void Initialize(size_t num_variables, size_t history_size, size_t history_subsampling);

		void Reset();
		void AddSample(const VectorReal& sample);

		MatrixReal GetHistory() const;
		MatrixReal GetHistory(std::vector<ptrdiff_t> variable_indices) const;
		MatrixReal GetEmpiricalCovariance() const;
		MatrixReal GetEmpiricalCorrelation() const;

		size_t GetSampleCount() const;
		Eigen::VectorXf GetHistorySample(ptrdiff_t i) const;

	private:
		Eigen::MatrixXf samples;		// Note samples are stored in columns
		ptrdiff_t sample_n;
		ptrdiff_t sample_n_s;
		ptrdiff_t sample_subsampling;

		friend class SampleHistoryClustering;
	};

}
