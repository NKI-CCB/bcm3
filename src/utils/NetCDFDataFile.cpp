#include "Utils.h"
#include "NetCDFDataFile.h"
#include <netcdf.h>

// Not very pretty but saves a lot of code copying
#define NC_HANDLE_ERROR(fn, msg, ...) { int result = fn; if (result != NC_NOERR) { LOGERROR(msg, __VA_ARGS__, result); return false; } }

namespace bcm3 {

template<> int NetCDFDataFile::CreateVariableSpec<int>					(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const { return nc_def_var(group, varname, NC_INT	, 1, dimidsp, varidp); }
template<> int NetCDFDataFile::CreateVariableSpec<unsigned int>			(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const { return nc_def_var(group, varname, NC_UINT	, 1, dimidsp, varidp); }
template<> int NetCDFDataFile::CreateVariableSpec<long long>			(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const { return nc_def_var(group, varname, NC_INT64	, 1, dimidsp, varidp); }
template<> int NetCDFDataFile::CreateVariableSpec<unsigned long long>	(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const { return nc_def_var(group, varname, NC_UINT64, 1, dimidsp, varidp); }
template<> int NetCDFDataFile::CreateVariableSpec<float>				(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const { return nc_def_var(group, varname, NC_FLOAT	, 1, dimidsp, varidp); }
template<> int NetCDFDataFile::CreateVariableSpec<double>				(int group, const char* varname, int ndims, const int* dimidsp, int* varidp) const { return nc_def_var(group, varname, NC_DOUBLE, 1, dimidsp, varidp); }

template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, int*					value) const { return nc_get_var1_int		(group, var, &dim1ix, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, unsigned int*		value) const { return nc_get_var1_uint		(group, var, &dim1ix, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, long*				value) const { return nc_get_var1_long		(group, var, &dim1ix, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, long long*			value) const { return nc_get_var1_longlong	(group, var, &dim1ix, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, unsigned long long*	value) const { return nc_get_var1_ulonglong	(group, var, &dim1ix, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, float*				value) const { 
	int res = nc_get_var1_float (group, var, &dim1ix, value);
	if (*value == NC_FILL_FLOAT) {
		*value = std::numeric_limits<float>::quiet_NaN();
	}
	return res;
}
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, double*				value) const {
	int res = nc_get_var1_double(group, var, &dim1ix, value);
	if (*value == NC_FILL_DOUBLE) {
		*value = std::numeric_limits<double>::quiet_NaN();
	}
	return res;
}

template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, int*					value) const { size_t dimids[2] = { dim1ix, dim2ix }; return nc_get_var1_int		(group, var, dimids, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, unsigned int*			value) const { size_t dimids[2] = { dim1ix, dim2ix }; return nc_get_var1_uint		(group, var, dimids, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, long*					value) const { size_t dimids[2] = { dim1ix, dim2ix }; return nc_get_var1_long		(group, var, dimids, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, long long*			value) const { size_t dimids[2] = { dim1ix, dim2ix }; return nc_get_var1_longlong	(group, var, dimids, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, unsigned long long*	value) const { size_t dimids[2] = { dim1ix, dim2ix }; return nc_get_var1_ulonglong	(group, var, dimids, value); }
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, float*				value) const {
	size_t dimids[2] = { dim1ix, dim2ix };
	int res = nc_get_var1_float (group, var, dimids, value);
	if (*value == NC_FILL_FLOAT) {
		*value = std::numeric_limits<float>::quiet_NaN();
	}
	return res;
}
template<> int NetCDFDataFile::GetValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, double*				value) const {
	size_t dimids[2] = { dim1ix, dim2ix };
	int res = nc_get_var1_double(group, var, dimids, value);
	if (*value == NC_FILL_DOUBLE) {
		*value = std::numeric_limits<double>::quiet_NaN();
	}
	return res;
}

template<> int NetCDFDataFile::GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<int>&					values) const { return nc_get_vara_int(group, var, &dim1ix, &count, values.data()); }
template<> int NetCDFDataFile::GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<unsigned int>&		values) const { return nc_get_vara_uint(group, var, &dim1ix, &count, values.data()); }
template<> int NetCDFDataFile::GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<long>&				values) const { return nc_get_vara_long(group, var, &dim1ix, &count, values.data()); }
template<> int NetCDFDataFile::GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<long long>&			values) const { return nc_get_vara_longlong(group, var, &dim1ix, &count, values.data()); }
template<> int NetCDFDataFile::GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<unsigned long long>&	values) const { return nc_get_vara_ulonglong(group, var, &dim1ix, &count, values.data()); }
template<> int NetCDFDataFile::GetValuesSpec(int group, int var, size_t dim1ix, size_t count, std::vector<std::string>&			values) const { 
	std::vector<char*> cvalues(count);
	int res = nc_get_vara_string(group, var, &dim1ix, &count, cvalues.data());
	if (res != 0) {
		return res;
	}
	for (size_t i = 0; i < count; i++) {
		values[i] = cvalues[i];
	}
	return 0;
}

template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, int					value) { return nc_put_var1_int			(group, var, &dim1ix, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, unsigned int			value) { return nc_put_var1_uint		(group, var, &dim1ix, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, long					value) { return nc_put_var1_long		(group, var, &dim1ix, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, long long			value) { return nc_put_var1_longlong	(group, var, &dim1ix, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, unsigned long long	value) { return nc_put_var1_ulonglong	(group, var, &dim1ix, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, float				value) { return nc_put_var1_float		(group, var, &dim1ix, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, double				value) { return nc_put_var1_double		(group, var, &dim1ix, &value); }

template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, int					value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_int		(group, var, dimids, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, unsigned int			value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_uint		(group, var, dimids, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, long					value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_long		(group, var, dimids, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, long long				value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_longlong	(group, var, dimids, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, unsigned long long	value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_ulonglong(group, var, dimids, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, float					value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_float	(group, var, dimids, &value); }
template<> int NetCDFDataFile::PutValueSpec(int group, int var, size_t dim1ix, size_t dim2ix, double				value) { size_t dimids[2] = { dim1ix, dim2ix }; return nc_put_var1_double	(group, var, dimids, &value); }

NetCDFDataFile::NetCDFDataFile()
	: fd(-1)
{
}

NetCDFDataFile::~NetCDFDataFile()
{
}

bool NetCDFDataFile::Open(const std::string& filename, bool write)
{
	int result;
	if (write) {
		result = nc_open(filename.c_str(), NC_WRITE, &fd);
	} else {
		result = nc_open(filename.c_str(), NC_NOWRITE, &fd);
	}
	if (result != NC_NOERR) {
		LOGERROR("Error opening netcdf file \"%s\"; status: %d", filename.c_str(), result);
		return false;
	} else {
		return true;
	}
}

bool NetCDFDataFile::Create(const std::string& filename)
{
	int result;
	result = nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &fd);
	if (result != NC_NOERR) {
		LOGERROR("Error creating netcdf file \"%s\"; status: %d", filename.c_str(), result);
		return false;
	} else {
		return true;
	}
}

void NetCDFDataFile::Close()
{
	int result = nc_close(fd);
	if (result != NC_NOERR) {
		LOGERROR("Error closing netcdf file; status: %d", result);
	}
}

void NetCDFDataFile::Sync()
{
	int result = nc_sync(fd);
	if (result != NC_NOERR) {
		LOGERROR("Error syncing netcdf file; status: %d", result);
	}
}

bool NetCDFDataFile::CreateGroup(const std::string& group_name)
{
	size_t slash = group_name.find('/');
	if (slash != std::string::npos && slash < group_name.size() + 1) {
		std::string subgroup = group_name.substr(0, slash);
		std::string rest = group_name.substr(slash + 1);

		int grpid, grpid2;
		int result = nc_inq_ncid(fd, subgroup.c_str(), &grpid);
		if (result == NC_NOERR) {
			int result = nc_inq_ncid(grpid, rest.c_str(), &grpid2);
			if (result == NC_NOERR) {
				return true;
			} else {
				NC_HANDLE_ERROR(nc_def_grp(grpid, rest.c_str(), &grpid2), "Error creating group \"%s\"; status: %d", group_name.c_str())
				return true;
			}
		} else {
			NC_HANDLE_ERROR(nc_def_grp(fd, subgroup.c_str(), &grpid), "Error creating group \"%s\"; status: %d", group_name.c_str())
			NC_HANDLE_ERROR(nc_def_grp(grpid, rest.c_str(), &grpid2), "Error creating group \"%s\"; status: %d", group_name.c_str())
			return true;
		}
	} else {
		int grpid;
		NC_HANDLE_ERROR(nc_def_grp(fd, group_name.c_str(), &grpid), "Error creating group \"%s\"; status: %d", group_name.c_str())
		return true;
	}
}

bool NetCDFDataFile::CreateDimension(const std::string& group_name, const std::string& dimname, const std::vector<unsigned int>& dimvalues)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }

	int dimid = -1, varid = -1;
	NC_HANDLE_ERROR(nc_def_dim(group, dimname.c_str(), dimvalues.size(), &dimid), "Error creating dimension \"%s\"; status: %d", dimname.c_str())
	NC_HANDLE_ERROR(nc_def_var(group, dimname.c_str(), NC_UINT, 1, &dimid, &varid), "Error creating dimension variable \"%s\"; status: %d", dimname.c_str())
	size_t start = 0;
	size_t count = dimvalues.size();
	NC_HANDLE_ERROR(nc_put_vara_uint(group, varid, &start, &count, dimvalues.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::CreateDimension(const std::string& group_name, const std::string& dimname, const std::vector<std::string>& dimvalues)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }

	int dimid = -1, varid = -1;
	NC_HANDLE_ERROR(nc_def_dim(group, dimname.c_str(), dimvalues.size(), &dimid), "Error creating dimension \"%s\"; status: %d", dimname.c_str())
	NC_HANDLE_ERROR(nc_def_var(group, dimname.c_str(), NC_STRING, 1, &dimid, &varid), "Error creating dimension variable \"%s\"; status: %d", dimname.c_str())
	for (size_t i = 0; i < dimvalues.size(); i++) {
		size_t count = 1;
		const char* p = dimvalues[i].c_str();
		NC_HANDLE_ERROR(nc_put_vara_string(group, varid, &i, &count, &p), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str())
	}
	return true;
}

bool NetCDFDataFile::CreateDimension(const std::string& group_name, const std::string& dimname, const VectorReal& dimvalues)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }

	int dimid = -1, varid = -1;
	NC_HANDLE_ERROR(nc_def_dim(group, dimname.c_str(), dimvalues.size(), &dimid), "Error creating dimension \"%s\"; status: %d", dimname.c_str())
	NC_HANDLE_ERROR(nc_def_var(group, dimname.c_str(), NC_DOUBLE, 1, &dimid, &varid), "Error creating dimension variable \"%s\"; status: %d", dimname.c_str())
	size_t start = 0;
	size_t count = dimvalues.size();
	NC_HANDLE_ERROR(nc_put_vara_double(group, varid, &start, &count, dimvalues.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int dim = GetDimension(group, dim1name);
	if (dim == -1) { return false; }

	int varid = -1;
	NC_HANDLE_ERROR(nc_def_var(group, variable_name.c_str(), NC_DOUBLE, 1, &dim, &varid), "Error creating variable \"%s\"; status: %d", variable_name.c_str())
	return true;
}

bool NetCDFDataFile::CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name, const std::string& dim2name)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }

	int dims[2];
	dims[0] = GetDimension(group, dim1name);
	if (dims[0] == -1) { return false; }
	dims[1] = GetDimension(group, dim2name);
	if (dims[1] == -1) { return false; }

	int varid = -1;
	NC_HANDLE_ERROR(nc_def_var(group, variable_name.c_str(), NC_DOUBLE, 2, dims, &varid), "Error creating variable \"%s\"; status: %d", variable_name.c_str())
	return true;
}

bool NetCDFDataFile::CreateVariable(const std::string& group_name, const std::string& variable_name, const std::string& dim1name, const std::string& dim2name, const std::string& dim3name)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }

	int dims[3];
	dims[0] = GetDimension(group, dim1name);
	if (dims[0] == -1) { return false; }
	dims[1] = GetDimension(group, dim2name);
	if (dims[1] == -1) { return false; }
	dims[2] = GetDimension(group, dim3name);
	if (dims[2] == -1) { return false; }

	int varid = -1;
	NC_HANDLE_ERROR(nc_def_var(group, variable_name.c_str(), NC_DOUBLE, 3, dims, &varid), "Error creating variable \"%s\"; status: %d", group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, Real* value) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t dimids[3] = { dim1ix, dim2ix, dim3ix };
	NC_HANDLE_ERROR(nc_get_var1_double(group, var, dimids, value), "Error getting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (*value == NC_FILL_FLOAT) {
		*value = std::numeric_limits<Real>::quiet_NaN();
	}
	return true;
}

bool NetCDFDataFile::GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t dim4ix, Real* value) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 4) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t dimids[4] = { dim1ix, dim2ix, dim3ix, dim4ix };
	NC_HANDLE_ERROR(nc_get_var1_double(group, var, dimids, value), "Error getting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (*value == NC_FILL_FLOAT) {
		*value = std::numeric_limits<Real>::quiet_NaN();
	}
	return true;
}

bool NetCDFDataFile::GetValue(const std::string& group_name, const std::string& variable_name, size_t dim1ix, std::string& value) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 1) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	char* string;
	NC_HANDLE_ERROR(nc_get_var1_string(group, var, &dim1ix, &string), "Error getting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (string != nullptr) {
		value = string;
#if _DEBUG
		// This will give problems when linking to "Release" netCDF dll as we can't mix memory allocation.. accept a memory leak instead?
#else
		free(string);
#endif
		return true;
	} else {
		return false;
	}
}

bool NetCDFDataFile::GetValue(const std::string& group_name, const std::string& variable_name, const std::vector<std::string>& ix, Real* value) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != ix.size()) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t* dimids = new size_t[ix.size()];
	for (size_t i = 0; i < ix.size(); i++) {
		if (!GetDimensionIx(group_name, variable_name, ix[i], &dimids[i])) {
			delete[] dimids;
			return false;
		}
	}

	NC_HANDLE_ERROR(nc_get_var1_double(group, var, dimids, value), "Error getting value for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (*value == NC_FILL_FLOAT) {
		*value = std::numeric_limits<Real>::quiet_NaN();
	}
	delete[] dimids;
	return true;
}

bool NetCDFDataFile::GetValues(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t count, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 1) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, &dim1ix, &count, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t count, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 2) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	size_t start[2] = { dim1ix, dim2ix };
	size_t ncount[2] = { count, 1 };
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, start, ncount, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t count, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 2) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	size_t start[2] = { dim1ix, dim2ix };
	size_t ncount[2] = { 1, count };
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, start, ncount, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t count, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	size_t start[3] = { dim1ix, dim2ix, dim3ix };
	size_t ncount[3] = { count, 1, 1 };
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, start, ncount, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t count, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	size_t start[3] = { dim1ix, dim2ix, dim3ix };
	size_t ncount[3] = { 1, count, 1 };
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, start, ncount, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetValuesDim3(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, size_t count, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to get values from variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	values.resize(count);
	size_t start[3] = { dim1ix, dim2ix, dim3ix };
	size_t ncount[3] = { 1, 1, count };
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, start, ncount, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetVector(const std::string& group_name, const std::string& variable_name, VectorReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 1) {
		LOGERROR("Trying to get vector from a variable that does not have one dimension, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	std::string dimname;
	bool result = true;
	result &= GetDimensionName(group_name, variable_name, 0, dimname);

	size_t dimsize;
	result &= GetDimensionSize(group_name, dimname, &dimsize);

	if (!result) {
		return false;
	}

	values.resize(dimsize);
	size_t dimix = 0;
	NC_HANDLE_ERROR(nc_get_vara_double(group, var, &dimix, &dimsize, values.data()), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())

	return true;
}

bool NetCDFDataFile::GetMatrix(const std::string& group_name, const std::string& variable_name, MatrixReal& values) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 2) {
		LOGERROR("Trying to get matrix from a variable that does not have two dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	std::string dimname1, dimname2;
	bool result = true;
	result &= GetDimensionName(group_name, variable_name, 0, dimname1);
	result &= GetDimensionName(group_name, variable_name, 1, dimname2);

	size_t dim1size, dim2size;
	result &= GetDimensionSize(group_name, dimname1, &dim1size);
	result &= GetDimensionSize(group_name, dimname2, &dim2size);

	if (!result) {
		return false;
	}

	values.resize(dim1size, dim2size);
	size_t dimix[2] = { 0, 0 };
	size_t count[2] = { 1, dim2size };
	for (size_t i = 0; i < dim1size; i++) {
		dimix[0] = i;
		NC_HANDLE_ERROR(nc_get_vara_double(group, var, dimix, count, values.data() + i * dim2size), "Error getting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	}

	return true;
}

bool NetCDFDataFile::PutValues(const std::string& group_name, const std::string& variable_name, size_t dim1ix, const VectorReal& values)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 1) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t start = dim1ix;
	size_t count = values.size();
	NC_HANDLE_ERROR(nc_put_vara_double(group, var, &start, &count, values.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::PutValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, const VectorReal& values)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 2) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t start[2] = { dim1ix, dim2ix };
	size_t count[2] = { (size_t)values.size(), 1 };
	NC_HANDLE_ERROR(nc_put_vara_double(group, var, start, count, values.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::PutValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, const VectorReal& values)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
		if (ndims != 2) {
			LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
			return false;
		}

	size_t start[2] = { dim1ix, dim2ix };
	size_t count[2] = { 1, (size_t)values.size() };
	NC_HANDLE_ERROR(nc_put_vara_double(group, var, start, count, values.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::PutValuesDim1(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, const VectorReal& values)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t start[3] = { dim1ix, dim2ix, dim3ix };
	size_t count[3] = { (size_t)values.size(), 1, 1 };
	NC_HANDLE_ERROR(nc_put_vara_double(group, var, start, count, values.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::PutValuesDim2(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, const VectorReal& values)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t start[3] = { dim1ix, dim2ix, dim3ix };
	size_t count[3] = { 1, (size_t)values.size(), 1 };
	NC_HANDLE_ERROR(nc_put_vara_double(group, var, start, count, values.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::PutValuesDim3(const std::string& group_name, const std::string& variable_name, size_t dim1ix, size_t dim2ix, size_t dim3ix, const VectorReal& values)
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (ndims != 3) {
		LOGERROR("Trying to write to variable with incorrect number of dimensions, for variable \"%s\" in group \"%s\"; status: %d", group_name.c_str(), variable_name.c_str());
		return false;
	}

	size_t start[3] = { dim1ix, dim2ix, dim3ix };
	size_t count[3] = { 1, 1, (size_t)values.size() };
	NC_HANDLE_ERROR(nc_put_vara_double(group, var, start, count, values.data()), "Error putting values for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::VariableExists(const std::string& group_name, const std::string& variable_name) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = -1;
	int result = nc_inq_varid(group, variable_name.c_str(), &var);
	return result == NC_NOERR;
}

bool NetCDFDataFile::GetDimensionCount(const std::string& group_name, const std::string& variable_name, size_t* count) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error retrieving number of dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	*count = ndims;
	return true;
}

bool NetCDFDataFile::GetDimensionName(const std::string& group_name, const std::string& variable_name, size_t dimix, std::string& dimname) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, variable_name);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error retrieving number of dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	if (dimix >= ndims) {
		LOGERROR("Dimension name requested for dimension index %d but the variable only has %d dimensions, variable \"%s\" in group \"%s\"; status: %d", dimix, ndims, variable_name.c_str(), group_name.c_str());
		return false;
	}

	int* dimids = new int[ndims];
	NC_HANDLE_ERROR(nc_inq_vardimid(group, var, dimids), "Error retrieving dimension ids for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str())
	char cdimname[NC_MAX_NAME];
	NC_HANDLE_ERROR(nc_inq_dimname(group, dimids[dimix], cdimname), "Error retrieving dimension name for dimension %d of variable \"%s\" in group \"%s\"; status: %d", dimix, variable_name.c_str(), group_name.c_str())
	dimname = cdimname;
	return true;
}

bool NetCDFDataFile::GetDimensionSize(const std::string& group_name, const std::string& dimname, size_t* size) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int dim = GetDimension(group, dimname);
	if (dim == -1) { return false; }

	NC_HANDLE_ERROR(nc_inq_dimlen(group, dim, size), "Error retrieving dimension size for dimension \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str())
	return true;
}

bool NetCDFDataFile::GetDimensionIx(const std::string& group_name, const std::string& dimname, const std::string& dimvalue, size_t* ix) const
{
	int group = GetGroup(group_name);
	if (group == -1) { return false; }
	int var = GetVariable(group, dimname);
	if (var == -1) { return false; }

	int ndims;
	NC_HANDLE_ERROR(nc_inq_varndims(group, var, &ndims), "Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str())
	if (ndims != 1) {
		LOGERROR("Dimensional variable has more than 1 dimension, for dimension \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str());
		return false;
	}

	int dimid;
	NC_HANDLE_ERROR(nc_inq_vardimid(group, var, &dimid), "Error retrieving dimension ids for variable \"%s\" in group \"%s\"; status: %d", dimvalue.c_str(), group_name.c_str())

	size_t dimsize;
	NC_HANDLE_ERROR(nc_inq_dimlen(group, dimid, &dimsize), "Error retrieving dimension size for dimension \"%s\" in group \"%s\"; status: %d", dimname.c_str(), group_name.c_str())

	for (size_t i = 0; i < dimsize; i++) {
		char* string;
		NC_HANDLE_ERROR(nc_get_var1_string(group, var, &i, &string), "Error retrieving value %d for dimension \"%s\" in group \"%s\"; status: %d", i, dimname.c_str(), group_name.c_str())
		if (string != nullptr) {
			if (strcmp(string, dimvalue.c_str()) == 0) {
				*ix = i;
#if _DEBUG
				// This will give problems when linking to "Release" netCDF dll as we can't mix memory allocation.. accept a memory leak instead?
#else
				free(string);
#endif
				return true;
			} else {
#if _DEBUG
				// This will give problems when linking to "Release" netCDF dll as we can't mix memory allocation.. accept a memory leak instead?
#else
				free(string);
#endif
			}
		} else {
			return false;
		}
	}
	return false;
}

int NetCDFDataFile::GetGroup(const std::string& group_name) const
{
#if 1
	std::vector<std::string> hierarchy;
	bcm3::tokenize(group_name, hierarchy, "/");

	int group = fd;
	for (size_t i = 0; i < hierarchy.size(); i++) {
		int nextgroup;
		int result = nc_inq_ncid(group, hierarchy[i].c_str(), &nextgroup);
		if (result == NC_NOERR) {
			group = nextgroup;
		} else {
			LOGERROR("Error retrieving group \"%s\"; status: %d", group_name.c_str(), result);
		}
	}
#else
	int group = -1;
	int result = nc_inq_ncid(fd, group_name.c_str(), &group);
	if (result != NC_NOERR) {
		LOGERROR("Error retrieving group \"%s\"; status: %d", group_name.c_str(), result);
	}
#endif

	return group;
}

int NetCDFDataFile::GetDimension(int group, const std::string& dimension_name) const
{
	int dim = -1;
	int result = nc_inq_dimid(group, dimension_name.c_str(), &dim);
	if (result != NC_NOERR) {
		LOGERROR("Error retrieving dimension \"%s\"; status: %d", dimension_name.c_str(), result);
	}
	return dim;
}

int NetCDFDataFile::GetVariable(int group, const std::string& variable_name) const
{
	int var = -1;
	int result = nc_inq_varid(group, variable_name.c_str(), &var);
	if (result != NC_NOERR) {
		LOGERROR("Error retrieving variable \"%s\"; status: %d", variable_name.c_str(), result);
	}
	return var;
}


bool NetCDFDataFile::GetNDims(const std::string& group_name, const std::string& variable_name, int group, int var, int* ndims) const
{
	int result = nc_inq_varndims(group, var, ndims);
	if (result != NC_NOERR) {
		LOGERROR("Error inquiring dimensions for variable \"%s\" in group \"%s\"; status: %d", variable_name.c_str(), group_name.c_str(), result);
		return false;
	}
	return true;
}

}
