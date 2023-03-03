#pragma once

namespace bcm3 {

Real dmvnormal(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = false);
Real pmvnormal(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = false);
VectorReal dmvnormal_array(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = false);
VectorReal pmvnormal_array(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = false);

}
