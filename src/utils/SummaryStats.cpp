#include "Utils.h"
#include "SummaryStats.h"

namespace bcm3 {

Real mean(const VectorReal& x)
{
	Real mu = 0.0;
	for (ptrdiff_t i = 0; i < x.size(); ++i) {
		mu += (x[i] - mu) / static_cast<Real>(i + 1);
	}
	return mu;
}

Real mean(const Eigen::VectorXf& x)
{
	Real mu = 0.0;
	for (ptrdiff_t i = 0; i < x.size(); ++i) {
		mu += (x[i] - mu) / static_cast<Real>(i + 1);
	}
	return mu;
}

Real var(const VectorReal& x)
{
	const Real mu = mean(x);
	return var(x, mu);
}

Real var(const Eigen::VectorXf& x)
{
	const Real mu = mean(x);
	return var(x, mu);
}

Real var(const VectorReal& x, const Real mu)
{
	Real sigmaSq = 0.0;
	for (ptrdiff_t i = 0; i < x.size(); ++i) {
		const Real d = x[i] - mu;
		sigmaSq += (d * d - sigmaSq) / static_cast<Real>(i + 1);
	}
	const Real n = static_cast<Real>(x.size());
	return sigmaSq * (n / (n - static_cast<Real>(1.0)));
}

Real var(const Eigen::VectorXf& x, const Real mu)
{
	Real sigmaSq = 0.0;
	for (ptrdiff_t i = 0; i < x.size(); ++i) {
		const Real d = x[i] - mu;
		sigmaSq += (d * d - sigmaSq) / static_cast<Real>(i + 1);
	}
	const Real n = static_cast<Real>(x.size());
	return sigmaSq * (n / (n - static_cast<Real>(1.0)));
}

Real acf(const VectorReal& x, const int lag)
{
	if (lag == 0) {
		return 1.0;
	}
	if (x.size() <= lag) {
		return std::numeric_limits<Real>::quiet_NaN();
	}

	const Real mu = mean(x);
	const Real sigmaSq = var(x, mu);

	Real r = 0.0;
	for (ptrdiff_t i = 0; i < x.size() - lag; i++) {
		const Real x1 = x[i] - mu;
		const Real x2 = x[i + lag] - mu;
		r += (x1 * x2 - r) / static_cast<Real>(i + 1);
	}

	return r / sigmaSq;
}

Real acf(const VectorReal& x, const int lag, const Real mu, const Real sigmaSq)
{
	if (lag == 0) {
		return 1.0;
	}
	if (x.size() <= lag) {
		return std::numeric_limits<Real>::quiet_NaN();
	}

	Real r = 0.0;
	for (ptrdiff_t i = 0; i < x.size() - lag; i++) {
		const Real x1 = x[i] - mu;
		const Real x2 = x[i + lag] - mu;
		r += (x1 * x2 - r) / static_cast<Real>(i + 1);
	}

	return r / sigmaSq;
}

VectorReal rowSum(const MatrixReal& x)
{
	VectorReal v(x.rows());
	for (ptrdiff_t i = 0; i < x.rows(); i++) {
		v(i) = x.row(i).sum();
	}
	return v;
}

VectorReal rowMean(const MatrixReal& x)
{
	VectorReal v(x.rows());
	for (ptrdiff_t i = 0; i < x.rows(); i++) {
		Real mu = 0.0;
		for (ptrdiff_t j = 0; j < x.cols(); ++j) {
			mu += (x(i,j) - mu) / static_cast<Real>(j + 1);
		}
		v(i) = mu;
	}
	return v;
}

VectorReal rowVar(const MatrixReal& x)
{
	VectorReal v(x.rows());
	for (ptrdiff_t i = 0; i < x.rows(); i++) {
		Real mu = 0.0;
		for (ptrdiff_t j = 0; j < x.cols(); ++j) {
			mu += (x(i, j) - mu) / static_cast<Real>(j + 1);
		}

		Real sigmaSq = 0.0;
		for (ptrdiff_t j = 0; j < x.cols(); ++j) {
			const Real d = x(i, j) - mu;
			sigmaSq += (d * d - sigmaSq) / static_cast<Real>(j + 1);
		}
		const Real n = static_cast<Real>(x.cols());

		v(i) = sigmaSq * (n / (n - static_cast<Real>(1.0)));
	}
	return v;
}

VectorReal rowSd(const MatrixReal& x)
{
	return rowVar(x).cwiseSqrt();
}

VectorReal colSum(const MatrixReal& x)
{
	VectorReal v(x.cols());
	for (ptrdiff_t i = 0; i < x.cols(); i++) {
		v(i) = x.col(i).sum();
	}
	return v;
}

VectorReal colMean(const MatrixReal& x)
{
	VectorReal v(x.cols());
	for (ptrdiff_t i = 0; i < x.cols(); i++) {
		Real mu = 0.0;
		for (ptrdiff_t j = 0; j < x.rows(); ++j) {
			mu += (x(j, i) - mu) / static_cast<Real>(j + 1);
		}
		v(i) = mu;
	}
	return v;
}

VectorReal colVar(const MatrixReal& x)
{
	VectorReal v(x.cols());
	for (ptrdiff_t i = 0; i < x.cols(); i++) {
		Real mu = 0.0;
		for (ptrdiff_t j = 0; j < x.rows(); ++j) {
			mu += (x(j, i) - mu) / static_cast<Real>(j + 1);
		}

		Real sigmaSq = 0.0;
		for (ptrdiff_t j = 0; j < x.rows(); ++j) {
			const Real d = x(j, i) - mu;
			sigmaSq += (d * d - sigmaSq) / static_cast<Real>(j + 1);
		}
		const Real n = static_cast<Real>(x.cols());

		v(i) = sigmaSq * (n / (n - static_cast<Real>(1.0)));
	}
	return v;
}

VectorReal colSd(const MatrixReal& x)
{
	return colVar(x).cwiseSqrt();
}

MatrixReal cov(const MatrixReal& x)
{
	MatrixReal c(x.cols(), x.cols());

	VectorReal empirical_mean = VectorReal::Zero(x.cols());
	VectorReal empirical_mean_nm1 = VectorReal::Zero(x.cols());
	MatrixReal empirical_correlation_calc = MatrixReal::Zero(x.cols(), x.cols());

	Real empirical_n = 0.0;
	for (ptrdiff_t si = 0; si < x.rows(); si++) {
		empirical_n += 1.0;

		Real invN = 1.0 / empirical_n;
		empirical_mean_nm1 = empirical_mean;
		for (ptrdiff_t i = 0; i < x.cols(); i++) {
			const Real v = (Real)x(si, i);
			const Real d = v - empirical_mean[i];
			empirical_mean[i] += d * invN;
		}

		if (empirical_n > 1) {
			const Real ratio = (empirical_n - 1) / (Real)empirical_n;
			for (ptrdiff_t i = 0; i < x.cols(); i++) {
				const Real dx = x(si, i) - empirical_mean_nm1[i];
				for (ptrdiff_t j = i; j < x.cols(); j++) {
					const Real dy = x(si, j) - empirical_mean_nm1[j];
					empirical_correlation_calc(i, j) += dx * dy * ratio;
				}
			}
		}
	}

	// Calculate empirical correlation
	for (ptrdiff_t i = 0; i < x.cols(); i++) {
		for (ptrdiff_t j = 0; j < x.cols(); j++) {
			if (j >= i) {
				c(i, j) = empirical_correlation_calc(i, j);
			} else {
				c(i, j) = empirical_correlation_calc(j, i);
			}
		}
	}
	c *= 1.0 / (empirical_n - 1.0);

	return c;
}

MatrixReal cor(const MatrixReal& x)
{
	MatrixReal c(x.cols(), x.cols());

	VectorReal empirical_mean = VectorReal::Zero(x.cols());
	VectorReal empirical_mean_nm1 = VectorReal::Zero(x.cols());
	MatrixReal empirical_correlation_calc = MatrixReal::Zero(x.cols(), x.cols());

	Real empirical_n = 0.0;
	for (ptrdiff_t si = 0; si < x.rows(); si++) {
		empirical_n += 1.0;

		Real invN = 1.0 / empirical_n;
		empirical_mean_nm1 = empirical_mean;
		for (ptrdiff_t i = 0; i < x.cols(); i++) {
			const Real v = (Real)x(si, i);
			const Real d = v - empirical_mean[i];
			empirical_mean[i] += d * invN;
		}

		if (empirical_n > 1) {
			const Real ratio = (empirical_n - 1) / (Real)empirical_n;
			for (ptrdiff_t i = 0; i < x.cols(); i++) {
				const Real dx = x(si, i) - empirical_mean_nm1[i];
				for (ptrdiff_t j = i; j < x.cols(); j++) {
					const Real dy = x(si, j) - empirical_mean_nm1[j];
					empirical_correlation_calc(i, j) += dx * dy * ratio;
				}
			}
		}
	}

	// Calculate empirical correlation
	for (ptrdiff_t i = 0; i < x.cols(); i++) {
		for (ptrdiff_t j = 0; j < x.cols(); j++) {
			if (j >= i) {
				c(i, j) = empirical_correlation_calc(i, j) / (sqrt(empirical_correlation_calc(i, i)) * sqrt(empirical_correlation_calc(j, j)));
			} else {
				c(i, j) = empirical_correlation_calc(j, i) / (sqrt(empirical_correlation_calc(i, i)) * sqrt(empirical_correlation_calc(j, j)));
			}
		}
	}

	return c;
}

}
