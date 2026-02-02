#include "Utils.h"
#include "hungarian.h"
#include "matching.h"

int main(int argc, char* argv[])
{
	MatrixReal test(3, 3);
	test.row(0) << 2.2,  5.4, 1.1;
	test.row(1) << -2.5, 3.3, 8.0;
	test.row(2) << 7.0,  0.5, 0.1;

	MatrixReal matching(3, 3);
	bool result = HungarianMatching(test, matching, HUNGARIAN_MATCH_MAX);

	std::cout << matching;

	std::vector<WeightedBipartiteEdge> edges;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			edges.push_back(WeightedBipartiteEdge(i, j, -test(i, j)));
		}
	}

	std::vector<int> matching2 = hungarianMinimumWeightPerfectMatching(3, edges, 9);

	return 0;
}
