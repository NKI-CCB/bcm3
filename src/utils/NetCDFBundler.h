#pragma once

#include "NetCDFDataFile.h"

namespace bcm3 {

class NetCDFBundler
{
public:
	NetCDFBundler();
	~NetCDFBundler();

	bool Open(const std::string& filename);
	bool Close();

	bool AddGroup(const std::string& group);
	bool AddVector(const std::string& group, const std::string& name, const VectorReal& vec);
	bool AddMatrix(const std::string& group, const std::string& name, const MatrixReal& mat);
	bool AddString(const std::string& group, const std::string& name, const std::string& str);

	template<typename T>
	bool AddVector(const std::string& group, const std::string& name, const std::vector<T> vec);

private:
	NetCDFDataFile ncfile;
};

template<typename T>
bool NetCDFBundler::AddVector(const std::string& group, const std::string& name, const std::vector<T> vec)
{
	bool result = true;

	std::vector<unsigned int> dimvals;
	dimvals.resize(vec.size());
	for (size_t i = 0; i < dimvals.size(); i++) {
		dimvals[i] = i + 1;
	}

	std::string dimname = name + "_dim";

	result &= ncfile.CreateDimension(group, dimname, dimvals);
	result &= ncfile.CreateVariable<T>(group, name, dimname);
	for (size_t i = 0; i < dimvals.size(); i++) {
		result &= ncfile.PutValue(group, name, i, vec[i]);
	}

	return result;
}

}
