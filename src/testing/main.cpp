#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "RNG.h"
#include <fstream>

int main(int argc, char* argv[])
{
	const int D = 3;
	Real values[6] = { 0.1338472934, -3.4556032604, 0.4084691508, -0.5560838449, -1.0494807311, -1.9352298262 };

	MatrixReal cholesky_L(D, D);
	int value_reference_i = 0;
	for (int i = 0; i < D; i++) {
		Real scale = values[i];
		cholesky_L(i, i) = exp(scale);
		for (int j = 0; j < i; j++) {
			cholesky_L(i, j) = values[3 + (value_reference_i++)];
		}
		for (int j = i + 1; j < D; j++) {
			cholesky_L(i, j) = 0.0;
		}
	}

	bcm3::RNG rng;

	std::ofstream file("test.tsv");

	for (int ni = 0; ni < 500; ni++) {
		Real sample[3] = { rng.GetReal(), rng.GetReal(), rng.GetReal() };

		// Normally distributed pseudorandom values
		VectorReal v(D);
		int ix = 0;
		for (int i = 0; i < D; i++) {
			v(i) = bcm3::QuantileNormal(sample[ix++], 0, 1);
		}

		VectorReal x = cholesky_L * v;

		file << x.transpose() << std::endl;
	}
	file.close();


	return 0;
}
