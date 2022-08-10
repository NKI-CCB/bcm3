#pragma once

namespace bcm3 {

Real mean(const VectorReal& x);
Real mean(const Eigen::VectorXf& x);
Real var(const VectorReal& x);
Real var(const Eigen::VectorXf& x);
Real var(const VectorReal& x, const Real mu);
Real var(const Eigen::VectorXf& x, const Real mu);
Real acf(const VectorReal& x, const int lag);
Real acf(const VectorReal& x, const int lag, const Real mu, const Real sigmaSq);

VectorReal rowSum(const MatrixReal& x);
VectorReal rowMean(const MatrixReal& x);
VectorReal rowVar(const MatrixReal& x);
VectorReal rowSd(const MatrixReal& x);

VectorReal colSum(const MatrixReal& x);
VectorReal colMean(const MatrixReal& x);
VectorReal colVar(const MatrixReal& x);
VectorReal colSd(const MatrixReal& x);

MatrixReal cov(const MatrixReal& x);
MatrixReal cor(const MatrixReal& x);

}
