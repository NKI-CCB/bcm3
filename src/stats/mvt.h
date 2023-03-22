#pragma once

namespace bcm3 {

Real dmvt(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = false);
Real pmvt(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = false);
VectorReal dmvt_array(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = false);
VectorReal pmvt_array(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, Real nu, bool return_log = false);

}
