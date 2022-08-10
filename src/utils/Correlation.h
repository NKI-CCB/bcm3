#pragma once

namespace bcm3 {
Real pearson_correlation(const VectorReal& x1, const VectorReal& x2);
Real pearson_correlation_weighted(const VectorReal& x1, const VectorReal& x2, const VectorReal& weights);
Real spearman_correlation(const VectorReal& x1, const VectorReal& x2, VectorReal& work1, VectorReal& work2);
Real spearman_correlation_weighted(const VectorReal& x1, const VectorReal& x2, const VectorReal& weights, VectorReal& work1, VectorReal& work2, VectorReal& work3);

// In-place, changes values and order of x1 and x2 (and the order of the weights for the weighted version)
Real spearman_correlation(VectorReal& x1, VectorReal& x2);
Real spearman_correlation_weighted(VectorReal& x1, VectorReal& x2, VectorReal& weights);

}
