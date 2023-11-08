#include "Utils.h"
#include "SampleHistory.h"
#include "SummaryStats.h"

namespace bcm3 {

	SampleHistory::SampleHistory()
		: sample_n(0)
		, sample_n_s(0)
		, sample_subsampling(1)
	{
	}

	SampleHistory::~SampleHistory()
	{
	}

	void SampleHistory::Initialize(size_t num_variables, size_t history_size, size_t history_subsampling)
	{
		samples.resize(num_variables, history_size);
		sample_n = 0;
		sample_n_s = 0;
		sample_subsampling = history_subsampling;
	}

	void SampleHistory::Reset()
	{
		sample_n = 0;
		sample_n_s = 0;
	}

	void SampleHistory::AddSample(const VectorReal& sample)
	{
		sample_n_s++;
		if (sample_n_s == sample_subsampling) {
			size_t sample_history_ix = sample_n;
			if (sample_history_ix >= samples.cols()) {
				// Wrap around and start filling from the start again
				sample_history_ix = sample_history_ix % samples.cols();
			}
			samples.col(sample_history_ix) = sample.cast<float>();
			sample_n++;
			sample_n_s = 0;
		}
	}

	MatrixReal SampleHistory::GetHistory() const
	{
		ptrdiff_t n_history_samples = GetSampleCount();
		MatrixReal sample_history_d(n_history_samples, samples.rows());
		if (n_history_samples > 0) {
			sample_history_d = samples.block(0, 0, samples.rows(), n_history_samples).transpose().cast<Real>();
		}
		return sample_history_d;
	}

	MatrixReal SampleHistory::GetHistory(std::vector<ptrdiff_t> variable_indices) const
	{
		ptrdiff_t n_history_samples = GetSampleCount();
		MatrixReal sample_history_d(n_history_samples, variable_indices.size());
		if (n_history_samples > 0) {
			for (ptrdiff_t i = 0; i < variable_indices.size(); i++) {
				sample_history_d.col(i) = samples.row(variable_indices[i]).segment(0, n_history_samples).cast<Real>();
			}
		}
		return sample_history_d;
	}

	MatrixReal SampleHistory::GetEmpiricalCovariance() const
	{
		ptrdiff_t n_history_samples = std::min((ptrdiff_t)samples.cols(), sample_n);
		MatrixReal sample_history_d = samples.block(0, 0, samples.rows(), n_history_samples).cast<Real>();
		return cov(sample_history_d.transpose());
	}

	MatrixReal SampleHistory::GetEmpiricalCorrelation() const
	{
		ptrdiff_t n_history_samples = std::min((ptrdiff_t)samples.cols(), sample_n);
		MatrixReal sample_history_d = samples.block(0, 0, samples.rows(), n_history_samples).cast<Real>();
		return cor(sample_history_d.transpose());
	}

	size_t SampleHistory::GetSampleCount() const
	{
		return std::min((ptrdiff_t)samples.cols(), sample_n);
	}

	Eigen::VectorXf SampleHistory::GetHistorySample(ptrdiff_t i) const
	{
		return samples.col(i);
	}

}
