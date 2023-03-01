#pragma once

namespace bcm3 {

Real dmvt(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = true);
Real pmvt(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = true);
VectorReal dmvt(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = true);
VectorReal pmvt(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = true);

}
