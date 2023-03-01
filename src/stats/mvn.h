#pragma once

namespace bcm3 {

Real dmvnormal(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = true);
Real pmvnormal(const VectorReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = true);
VectorReal dmvnormal(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = true);
VectorReal pmvnormal(const MatrixReal& x, const VectorReal& mu, const MatrixReal& sigma, bool return_log = true);

}
