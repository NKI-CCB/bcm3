#include "Utils.h"
#include "Correlation.h"

#include <numeric>

namespace bcm3 {

void sift(VectorReal& x, size_t start, size_t end)
{
	Real v = x[start];

	while (start <= end/2) {
		size_t j = 2 * start;
		if (j < end && x[j] < x[j+1]) {
			j++;
		}
		if (!(v < x[j])) {
			break;
		}
		x[start] = x[j];
		start = j;
	}

	x[start] = v;
}

void sift_reverse(VectorReal& x, size_t start, size_t end)
{
	Real v = x[start];

	while (start <= end / 2) {
		size_t j = 2 * start;
		if (j < end && x[j] > x[j + 1]) {
			j++;
		}
		if (!(v > x[j])) {
			break;
		}
		x[start] = x[j];
		start = j;
	}

	x[start] = v;
}

void sort(VectorReal& x, bool descend)
{
	if (x.size() == 0) {
		return;
	}

	size_t end = x.size() - 1;
	size_t start = end / 2 + 1;

	if (descend) {
		do {
			start--;
			sift_reverse(x, start, end);
		} while (start > 0);

		while (end > 0) {
			std::swap(x[0], x[end]);
			end--;
			sift_reverse(x, 0, end);
		}
	} else {
		do {
			start--;
			sift(x, start, end);
		} while (start > 0);

		while (end > 0) {
			std::swap(x[0], x[end]);
			end--;
			sift(x, 0, end);
		}
	}
}

void sift(VectorReal& x, VectorReal& aux, size_t start, size_t end)
{
	Real v = x[start];
	Real vaux = aux[start];

	while (start <= end/2) {
		size_t j = 2 * start;
		if (j < end && x[j] < x[j+1]) {
			j++;
		}
		if (!(v < x[j])) {
			break;
		}
		x[start] = x[j];
		aux[start] = aux[j];
		start = j;
	}

	x[start] = v;
	aux[start] = vaux;
}

void sort(VectorReal& x, VectorReal& aux)
{
	if (x.size() == 0) {
		return;
	}

	size_t end = x.size() - 1;
	size_t start = end / 2 + 1;

	do {
		start--;
		sift(x, aux, start, end);
	} while (start > 0);

	while (end > 0) {
		std::swap(x[0], x[end]);
		std::swap(aux[0], aux[end]);
		end--;
		sift(x, aux, 0, end);
	}
}

void sift(VectorReal& x, VectorReal& aux1, VectorReal& aux2, size_t start, size_t end)
{
	Real v = x[start];
	Real vaux1 = aux1[start];
	Real vaux2 = aux2[start];

	while (start <= end/2) {
		size_t j = 2 * start;
		if (j < end && x[j] < x[j+1]) {
			j++;
		}
		if (!(v < x[j])) {
			break;
		}
		x[start] = x[j];
		aux1[start] = aux1[j];
		aux2[start] = aux2[j];
		start = j;
	}

	x[start] = v;
	aux1[start] = vaux1;
	aux2[start] = vaux2;
}

void sort(VectorReal& x, VectorReal& aux1, VectorReal& aux2)
{
	if (x.size() == 0) {
		return;
	}

	size_t end = x.size() - 1;
	size_t start = end / 2 + 1;

	do {
		start--;
		sift(x, aux1, aux2, start, end);
	} while (start > 0);

	while (end > 0) {
		std::swap(x[0], x[end]);
		std::swap(aux1[0], aux1[end]);
		std::swap(aux2[0], aux2[end]);
		end--;
		sift(x, aux1, aux2, 0, end);
	}
}

std::vector<ptrdiff_t> rank(const VectorReal& x)
{
	std::vector<ptrdiff_t> result(x.size());
	std::vector<ptrdiff_t> ix(x.size());
	std::iota(ix.begin(), ix.end(), 0);
	std::sort(ix.begin(), ix.end(), [&x](ptrdiff_t i1, ptrdiff_t i2) {return x(i1) < x(i2); });
	for (ptrdiff_t i = 0; i < x.size(); ++i) {
		result[ix[i]] = i+1;
	}
	return result;
}

ptrdiff_t minindex(const VectorReal& x)
{
	ptrdiff_t minix = -1;
	Real minval = std::numeric_limits<Real>::infinity();
	for (ptrdiff_t i = 0; i < x.size(); i++) {
		if (x(i) < minval) {
			minval = x(i);
			minix = i;
		}
	}
	return minix;
}

VectorReal cumsum(VectorReal& x)
{
	VectorReal res(x.size());
	Real sum = 0.0;
	for (int i = 0; i < x.size(); i++) {
		sum += x(i);
		res(i) = sum;
	}
	return res;
}

bool ParseVectorFromString(const std::string& str, VectorReal& v)
{
	std::vector<std::string> vec;
	bcm3::tokenize(str, vec, ";");
	v.resize(vec.size());
	for (size_t i = 0; i < vec.size(); i++) {
		try {
			Real value = boost::lexical_cast<Real>(vec[i]);
			v(i) = value;
		} catch (const boost::bad_lexical_cast& e) {
			LOGERROR("Could not cast value \"%s\" to a constant real value: %s", vec[i].c_str(), e.what());
			return false;
		}
	}
	return true;
}

std::string VectorToString(const VectorReal& x)
{
	if (x.size() > 0) {
		std::string s = std::to_string(x(0));
		for (size_t i = 1; i < x.size(); i++) {
			s += ",";
			s += std::to_string(x(i));
		}
		return s;
	} else {
		return std::string();
	}
}

std::string MatrixToString(const MatrixReal& x)
{
	if (x.cols() > 0 && x.rows() > 0) {
		std::string s;
		for (size_t j = 0; j < x.rows(); j++) {
			if (j > 0) {
				s += "\n";
			}
			s += std::to_string(x(j,0));
			for (size_t i = 1; i < x.cols(); i++) {
				s += ",";
				s += std::to_string(x(j,i));
			}
		}
		return s;
	} else {
		return std::string();
	}
}

}
