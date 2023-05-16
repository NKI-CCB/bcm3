#pragma once

namespace bcm3 {

class NetCDFDataFile
{
public:
	NetCDFDataFile();
	~NetCDFDataFile();

	bool Open(const std::string& filename, bool write);
	bool Create(const std::string& filename);
	void Close();
	void Sync();

	bool CreateGroup(const std::string& group_name);

	bool CreateDimension(const std::string& group_name, const std::string& dimname, const std::vector<unsigned int>& dimvalues);
	bool CreateDimension(const std::string& group_name, const std::string& dimname, const std::vector<std::string>& dimvalues);
	bool CreateDimension(const std::string& group_name, const std::string& dimname, const VectorReal& dimvalues);

	bool CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name);
	bool CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name, const std::string& dim2name);
	bool CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name, const std::string& dim2name, const std::string& dim3name);

	template<typename T>
	bool CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name);

	template<typename T>
	bool GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, T* value) const;
	template<typename T>
	bool GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, T* value) const;
	bool GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, Real* value) const;
	bool GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t dim4ix, Real* value) const;
	bool GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, std::string& value) const;
	bool GetValue(const std::string& group_name, const std::string& variable_name, const std::vector<std::string>& ix, Real* value) const;

	template<typename T>
	bool GetValues(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t count, std::vector<T>& values) const;
	bool GetValues(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t count, VectorReal& values) const;
	bool GetValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t count, VectorReal& values) const;
	bool GetValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t count, VectorReal& values) const;
	bool GetValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t count, VectorReal& values) const;
	bool GetValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t count, VectorReal& values) const;
	bool GetValuesDim3(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t count, VectorReal& values) const;

	bool GetVector(const std::string& group_name, const std::string& variable_name, VectorReal& values) const;
	bool GetMatrix(const std::string& group_name, const std::string& variable_name, MatrixReal& values) const;

	template<typename T>
	bool PutValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, T value);
	template<typename T>
	bool PutValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, T value);
	bool PutValues(const std::string& group_name, const std::string& variable_name, size_t dim1ix, const VectorReal& values);
	bool PutValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, const VectorReal& values);
	bool PutValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, const VectorReal& values);
	bool PutValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, const VectorReal& values);
	bool PutValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, const VectorReal& values);
	bool PutValuesDim3(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, const VectorReal& values);

	bool VariableExists(const std::string& group_name, const std::string& variable_name) const;
	bool GetDimensionCount(const std::string& group_name, const std::string& variable_name, size_t* count) const;
	bool GetDimensionName(const std::string& group_name, const std::string& variable_name, size_t dimix, std::string& dimname) const;
	bool GetDimensionSize(const std::string& group_name, const std::string& dimname, size_t* size) const;
	bool GetDimensionIx(const std::string& group_name, const std::string& dimname, const std::string& dimvalue, size_t* ix) const;

private:
	int GetGroup(const std::string& group_name) const;
	int GetDimension(int group, const std::string& dimension_name) const;
	int GetVariable(int group, const std::string& variable_name) const;
	bool GetNDims(const std::string& group_name, const std::string& variable_name, int group, int var, int* ndims) const;

	template<typename T> int CreateVariableSpec(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const;
	template<typename T> int GetValueSpec(int group, int var, size_t dim1ix, T* value) const;
	template<typename T> int GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, T* value) const;
	template<typename T> int GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<T>& value) const;
	template<typename T> int PutValueSpec(int group, int var, size_t dim1ix, T value);
	template<typename T> int PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, T value);

	int fd;
};

template<typename T>
bool NetCDFDataFile::CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int dim = GetDimension(group, dim1name);
	if (dim == -1) { return false; }

	int varid = -1;
	int result = CreateVariableSpec<T>(group, variable_name.c_str(), 1, &dim, &varid);
	if (result != 0) {
		LOGERROR("Error creating variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str(), result);
		return false;
	}
	return true;
}

template<typename T>
bool NetCDFDataFile::GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, T* value) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	if (!GetNDims(group_name, variable_name, group, var, &ndims)) {
		return false;
	}
	if (ndims != 1) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	int res = GetValueSpec<T>(group, var, dim1ix, value);
	if (res != 0) {
		LOGERROR("Error getting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str());
		return false;
	}
	return true;
}

template<typename T>
bool NetCDFDataFile::GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, T* value) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	if (!GetNDims(group_name, variable_name, group, var, &ndims)) {
		return false;
	}
	if (ndims != 2) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	int res = GetValueSpec<T>(group, var, dim1ix, dim2ix, value);
	if (res != 0) {
		LOGERROR("Error getting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str());
		return false;
	}
	return true;
}

template<typename T>
bool NetCDFDataFile::GetValues(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t count, std::vector<T>& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	if (!GetNDims(group_name, variable_name, group, var, &ndims)) {
		return false;
	}
	if (ndims != 1) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	int res = GetValuesSpec<T>(group, var, dim1ix, count, values);
	if (res != 0) {
		LOGERROR("Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str());
		return false;
	}
	return true;
}

template<typename T>
bool NetCDFDataFile::PutValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, T value)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	if (!GetNDims(group_name, variable_name, group, var, &ndims)) {
		return false;
	}
	if (ndims != 1) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	int res = PutValueSpec<T>(group, var, dim1ix, value);
	if (res != 0) {
		LOGERROR("Error putting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str());
		return false;
	}
	return true;
}

template<typename T>
bool NetCDFDataFile::PutValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, T value)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	if (!GetNDims(group_name, variable_name, group, var, &ndims)) {
		return false;
	}
	if (ndims != 2) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	int res = PutValueSpec<T>(group, var, dim1ix, dim2ix, value);
	if (res != 0) {
		LOGERROR("Error putting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str());
		return false;
	}
	return true;
}

}
