#include "Utils.h"
#include "Correlation.h"
#include "VectorUtils.h"

namespace bcm3 {
	
Real pearson_correlation(const VectorReal& x1, const VectorReal& x2)
{
	ASSERT(x1.size() == x2.size());

	Real sumsq1 = 0.0;
	Real sumsq2 = 0.0;
	Real sumc = 0.0;

	Real d1, d2;
	Real mu1 = 0.0, mu2 = 0.0;

	const size_t n = (size_t)x1.size();
	for (size_t i = 0; i < n; i++) {
		Real f = (Real)i / ((Real)i + (Real)1.0);
		d1 = x1[i] - mu1;
		d2 = x2[i] - mu2;
		sumsq1 += d1 * d1 * f;
		sumsq2 += d2 * d2 * f;
		sumc += d1 * d2 * f;
		f = (Real)1.0 / ((Real)i + (Real)1.0);
		mu1 += d1 * f;
		mu2 += d2 * f;
	}

	Real r = sumc / (sqrt(sumsq1) * sqrt(sumsq2));
	return r;
}

Real pearson_correlation_weighted(const VectorReal& x1, const VectorReal& x2, const VectorReal& weights)
{
	ASSERT(x1.size() == x2.size());
	ASSERT(x1.size() == weights.size());

	Real sumsq1 = 0.0;
	Real sumsq2 = 0.0;
	Real sumc = 0.0;

	Real d1, d2;
	Real mu1 = 0.0, mu2 = 0.0;

	Real acc_weight = 0;
	const size_t n = (size_t)x1.size();
	for (size_t i = 0; i < n; i++) {
		const Real w = weights[i];
		if (w <= 0.0) {
			continue;
		}
		Real f = w * acc_weight / (acc_weight + w);
		d1 = x1[i] - mu1;
		d2 = x2[i] - mu2;
		sumsq1 += d1 * d1 * f;
		sumsq2 += d2 * d2 * f;
		sumc += d1 * d2 * f;
		f = w / (acc_weight + w);
		mu1 += d1 * f;
		mu2 += d2 * f;
		acc_weight += w;
	}

	Real r = sumc / (sqrt(sumsq1) * sqrt(sumsq2));
	return r;
}

void compute_rank(VectorReal& x)
{
	ASSERT(x.size() >= 1);
	size_t i = 0;
	const size_t n = (size_t)x.size();
	while (i < n - 1) {
		if (x[i] == x[i+1]) {
			// Tie
			size_t j = i + 2;
			while (j < n && x[i] == x[j]) {
				j++;
			}

			Real rank = (Real)0.0;
			for (size_t k = i; k < j; k++) {
				rank += k + (Real)1.0;
			}
			rank /= (Real)(j - i);

			for (size_t k = i; k < j; k++) {
				x[k] = rank;
			}
			i = j;
		} else {
			x[i] = (Real)i + (Real)1.0;
			i++;
		}
	}
	if (i == x.size() - 1) {
		x[i] = (Real)x.size();
	}
}

Real spearman_correlation(const VectorReal& x1, const VectorReal& x2, VectorReal& work1, VectorReal& work2)
{
	ASSERT(x1.size() == x2.size());

	work1 = x1;
	work2 = x2;

	sort(work1, work2);
	compute_rank(work1);
	sort(work2, work1);
	compute_rank(work2);

	return pearson_correlation(work1, work2);
}

Real spearman_correlation_weighted(const VectorReal& x1, const VectorReal& x2, const VectorReal& weights, VectorReal& work1, VectorReal& work2, VectorReal& work3)
{
	ASSERT(x1.size() == x2.size());

	work1 = x1;
	work2 = x2;
	work3 = weights;

	sort(work1, work2, work3);
	compute_rank(work1);
	sort(work2, work1, work3);
	compute_rank(work2);

	return pearson_correlation_weighted(work1, work2, work3);
}

Real spearman_correlation(VectorReal& x1, VectorReal& x2)
{
	ASSERT(x1.size() == x2.size());

	sort(x1, x2);
	compute_rank(x1);
	sort(x2, x1);
	compute_rank(x2);

	return pearson_correlation(x1, x2);
}

Real spearman_correlation_weighted(VectorReal& x1, VectorReal& x2, VectorReal& weights)
{
	ASSERT(x1.size() == x2.size());

	sort(x1, x2, weights);
	compute_rank(x1);
	sort(x2, x1, weights);
	compute_rank(x2);

	return pearson_correlation_weighted(x1, x2, weights);
}

void linear_regress_columns(const MatrixReal::ConstColXpr& x, const MatrixReal::ConstColXpr& y, Real& offset, Real& scale)
{
	ASSERT(x.size() == y.size());

	Real mu_x = 0.0;
	Real mu_y = 0.0;
	Real xvar_calc = 0.0;
	Real cov_calc = 0.0;

	Real empirical_n = 0.0;
	for (ptrdiff_t i = 0; i < x.size(); i++) {
		if (std::isnan(x(i)) || std::isnan(y(i))) {
			continue;
		}

		empirical_n += 1.0;

		const Real invN = 1.0 / empirical_n;
		Real mu_x_nm1 = mu_x;
		Real mu_y_nm1 = mu_y;

		mu_x += (x(i) - mu_x) * invN;
		mu_y += (y(i) - mu_y) * invN;

		if (empirical_n > 1) {
			const Real ratio = (empirical_n - 1) / (Real)empirical_n;

			const Real dx = x(i) - mu_x_nm1;
			const Real dy = y(i) - mu_y_nm1;

			xvar_calc += dx * dx * ratio;
			cov_calc += dx * dy * ratio;
		}
	}

	if (empirical_n >= 2) {
		// No need to calculate the actual covariance; the scaled value suffices since we only need the ratio
		scale = cov_calc / xvar_calc;
		offset = mu_y - mu_x * scale;
	} else {
		// Don't set scale and offset
	}
}

}
