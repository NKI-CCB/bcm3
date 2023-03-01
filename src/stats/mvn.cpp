#include "Utils.h"
#include "checks.h"
#include "mvn.h"

namespace bcm3 {

Real dmvnormal(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log)
{
	ASSERT(x.size() == mu.size());
	ASSERT(mu.size() == sigma.rows());
	ASSERT(mu.size() == sigma.cols());
	ASSERT(is_positive_semi_definite(sigma));

	const int p = mu.size();

	Eigen::LLT<MatrixReal> llt;
	llt.compute(sigma);
	ASSERT(llt.info() == Eigen::Success);

	Real det = 0.0;
	for (int i = 0; i < p; i++) {
		det += log(llt.matrixL()(i, i));
	}
	Real logC = -det - 0.5 * p * log(2.0 * M_PI);

	VectorReal v = x - mu;
	llt.matrixL().solveInPlace(v);
	Real logp = logC - 0.5 * v.dot(v);

	if (return_log) {
		return logp;
	} else {
		return exp(logp);
	}
}

Real pmvnormal(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log)
{
	ASSERT(x.size() == mu.size());
	ASSERT(mu.size() == sigma.rows());
	ASSERT(mu.size() == sigma.cols());
	ASSERT(is_positive_semi_definite(sigma));

	ASSERT(false);
	return 0.0;
}

VectorReal dmvnormal(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log)
{
	ASSERT(x.cols() == mu.size());
	ASSERT(mu.size() == sigma.rows());
	ASSERT(mu.size() == sigma.cols());
	ASSERT(is_positive_semi_definite(sigma));

	const int p = mu.size();

	Eigen::LLT<MatrixReal> llt;
	llt.compute(sigma);
	ASSERT(llt.info() == Eigen::Success);

	Real det = 0.0;
	for (int i = 0; i < p; i++) {
		det += log(llt.matrixL()(i, i));
	}
	Real logC = -det - 0.5 * p * log(2.0 * M_PI);

	VectorReal logp(x.rows());
	for (int i = 0; i < x.rows(); i++) {
		VectorReal v = x.row(i) - mu;
		llt.matrixL().solveInPlace(v);
		logp(i) = logC - 0.5 * v.dot(v);
	}

	if (return_log) {
		return logp;
	} else {
		return logp.array().exp();
	}
}

VectorReal pmvnormal(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log)
{
	ASSERT(x.size() == mu.size());
	ASSERT(mu.size() == sigma.rows());
	ASSERT(mu.size() == sigma.cols());
	ASSERT(is_positive_semi_definite(sigma));

	ASSERT(false);
	VectorReal logp(x.rows());
	return logp;
}

}
