#include "Utils.h"
#include "HungarianMatching.h"

namespace bcm3 {

std::optional<Eigen::VectorXi> RectangularHungarian(const Eigen::MatrixXd& cost)
{
	if (cost.rows() > cost.cols()) {
		LOGERROR("Hungarian matching cost matrix should have at least as many columns as rows");
		return std::nullopt;
	}
	if (cost.rows() > std::numeric_limits<int>::max() || cost.cols() > std::numeric_limits<int>::max()) {
		LOGERROR("Hungarian matching cost matrix too large for VectorXi indices");
		return std::nullopt;
	}

	const int n = static_cast<int>(cost.rows());
	const int m = static_cast<int>(cost.cols());

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (std::isnan(cost(i, j)) || cost(i, j) == -std::numeric_limits<Real>::infinity()) {
				LOGERROR("Costs must be finite or +infinity");
				return std::nullopt;
			}
		}
	}

	std::vector<Real> u(n + 1), v(m + 1);
	std::vector<int> p(m + 1), way(m + 1);

	for (int i = 1; i <= n; ++i) {
		p[0] = i;
		std::vector<Real> slack(m + 1, std::numeric_limits<Real>::infinity());
		std::vector<bool> used(m + 1, false);
		int j0 = 0;

		do {
			used[j0] = true;
			const int i0 = p[j0];
			Real delta = std::numeric_limits<Real>::infinity();
			int j1 = -1;

			for (int j = 1; j <= m; ++j) {
				if (used[j]) {
					continue;
				}

				const double c = cost(i0 - 1, j - 1);
				if (c != std::numeric_limits<Real>::infinity()) {
					const Real reduced = static_cast<Real>(c) - u[i0] - v[j];
					if (!std::isfinite(reduced)) {
						LOGWARNING("Potential arithmetic overflow");
						return std::nullopt;
					}

					if (reduced < slack[j]) {
						slack[j] = reduced;
						way[j] = j0;
					}
				}

				if (slack[j] < delta) {
					delta = slack[j];
					j1 = j;
				}
			}

			if (j1 == -1) {
				return std::nullopt;
			}

			for (int j = 0; j <= m; ++j) {
				if (used[j]) {
					u[p[j]] += delta;
					v[j] -= delta;
				} else if (slack[j] != std::numeric_limits<Real>::infinity()) {
					slack[j] -= delta;
				}
			}
			j0 = j1;
		} while (p[j0] != 0);

		do {
			const int previous = way[j0];
			p[j0] = p[previous];
			j0 = previous;
		} while (j0 != 0);
	}

	Eigen::VectorXi match(n);
	for (int j = 1; j <= m; ++j) {
		if (p[j] != 0) {
			match(p[j] - 1) = j - 1;
		}
	}
	return match;
}

}
