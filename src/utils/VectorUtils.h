#pragma once

namespace bcm3 {

void sort(VectorReal& x, bool descend = false);
void sort(VectorReal& x, VectorReal& aux);
void sort(VectorReal& x, VectorReal& aux1, VectorReal& aux2);
std::vector<ptrdiff_t> order(const VectorReal& x);
ptrdiff_t minindex(const VectorReal& x);

VectorReal cumsum(VectorReal& x);

inline void remove_element(VectorReal& x, size_t ix)
{
	size_t new_n = x.size() - 1;
	x.segment(ix, new_n - ix) = x.segment(ix + 1, new_n - ix);
	x.conservativeResize(new_n);
}

template<typename T>
bool ParseSquareMatrixFromString(const std::string& str, Eigen::Matrix<T, -1, -1>& matrix, size_t N = 0)
{
	std::vector<std::string> vec;
	bcm3::tokenize(str, vec, ";");
	if (vec.size() != N * N) {
		LOGERROR("Dimension of matrix does not match the required size");
		return false;
	}
	matrix.resize(N, N);
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			T value = boost::lexical_cast<T>(vec[j*N+i]);
			matrix(i, j) = value;
		}
	}
	return true;
}

template<typename T>
bool ParseMatrixFromString(const std::string& str, Eigen::Matrix<T, -1, -1>& matrix)
{
	std::vector<std::string> rows;
	bcm3::tokenize(str, rows, ";");

	std::vector<std::string> row;
	bcm3::tokenize(rows[0], row, ",");

	matrix.resize(rows.size(), row.size());

	try {
		for (size_t i = 0; i < rows.size(); i++) {
			bcm3::tokenize(rows[i], row, ",");
			if (row.size() != matrix.cols()) {
				LOGERROR("Inconsistent matrix");
				return false;
			}
			for (size_t j = 0; j < row.size(); j++) {
				T value = boost::lexical_cast<T>(row[j]);
				matrix(i, j) = value;
			}
		}
	} catch (const boost::bad_lexical_cast& e) {
		LOGERROR("Could not cast string to a Real value: %s", e.what());
		return false;
	}
	return true;
}

bool ParseVectorFromString(const std::string& str, VectorReal& v);

std::string VectorToString(const VectorReal& x);
std::string MatrixToString(const MatrixReal& x);

}
