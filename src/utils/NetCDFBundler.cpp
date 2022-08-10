#include "NetCDFBundler.h"
#include <boost/filesystem.hpp>

namespace bcm3 {

NetCDFBundler::NetCDFBundler()
{
}

NetCDFBundler::~NetCDFBundler()
{
}

bool NetCDFBundler::Open(const std::string& filename)
{
	if (boost::filesystem::exists(filename)) {
		return ncfile.Open(filename, true);
	} else {
		return ncfile.Create(filename);
	}
}

bool NetCDFBundler::Close()
{
	ncfile.Close();
	return true;
}

bool NetCDFBundler::AddGroup(const std::string& group)
{
	return ncfile.CreateGroup(group);
}

bool NetCDFBundler::AddVector(const std::string& group, const std::string& name, const VectorReal& vec)
{
	bool result = true;

	std::vector<unsigned int> dimvals;
	dimvals.resize(vec.size());
	for (size_t i = 0; i < dimvals.size(); i++) {
		dimvals[i] = i + 1;
	}

	std::string dimname = name + "_dim";

	result &= ncfile.CreateDimension(group, dimname, dimvals);
	result &= ncfile.CreateVariable(group, name, dimname);
	result &= ncfile.PutValues(group, name, 0, vec);

	return result;
}

bool NetCDFBundler::AddMatrix(const std::string& group, const std::string& name, const MatrixReal& mat)
{
	bool result = true;

	std::string dimname1 = name + "_dim1";
	std::string dimname2 = name + "_dim2";

	std::vector<unsigned int> dimvals;
	dimvals.resize(mat.rows());
	for (size_t i = 0; i < dimvals.size(); i++) {
		dimvals[i] = i + 1;
	}
	result &= ncfile.CreateDimension(group, dimname1, dimvals);

	dimvals.resize(mat.cols());
	for (size_t i = 0; i < dimvals.size(); i++) {
		dimvals[i] = i + 1;
	}
	result &= ncfile.CreateDimension(group, dimname2, dimvals);

	result &= ncfile.CreateVariable(group, name, dimname1, dimname2);
	for (ptrdiff_t i = 0; i < mat.rows(); i++) {
		result &= ncfile.PutValuesDim2(group, name, i, 0, mat.row(i));
	}

	return result;
}

bool NetCDFBundler::AddString(const std::string& group, const std::string& name, const std::string& str)
{
	return true;
}

}