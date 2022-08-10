#include "Utils.h"
#include "DataReference.h"

namespace bcm3 {

std::shared_ptr<DataReference> DataReference::Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::string dim1, std::string ix1)
{
	std::vector<std::string> dims(1);
	dims[0] = dim1;
	std::vector<std::string> ixs(1);
	ixs[0] = ix1;
	return Create(data_file, group, variable_name, dims, ixs);
}

std::shared_ptr<DataReference> DataReference::Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::string dim1, std::string ix1, std::string dim2, std::string ix2)
{
	std::vector<std::string> dims(2);
	dims[0] = dim1;
	dims[1] = dim2;
	std::vector<std::string> ixs(2);
	ixs[0] = ix1;
	ixs[1] = ix2;
	return Create(data_file, group, variable_name, dims, ixs);
}

std::shared_ptr<DataReference> DataReference::Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::string dim1, std::string ix1, std::string dim2, std::string ix2, std::string dim3, std::string ix3)
{
	std::vector<std::string> dims(3);
	dims[0] = dim1;
	dims[1] = dim2;
	dims[2] = dim3;
	std::vector<std::string> ixs(3);
	ixs[0] = ix1;
	ixs[1] = ix2;
	ixs[2] = ix3;
	return Create(data_file, group, variable_name, dims, ixs);
}

std::shared_ptr<DataReference> DataReference::Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::vector<std::string> dimensions, std::vector<std::string> indices)
{
	std::shared_ptr<DataReference> ref;

	if (dimensions.size() != indices.size()) {
		LOGERROR("Inconsistent dimensions/indices for data reference to %s/%s; %d/%d", group.c_str(), variable_name.c_str(), dimensions.size(), indices.size());
		return ref;
	}

	size_t dimcount;
	if (!data_file.GetDimensionCount(group, variable_name, &dimcount)) {
		return ref;
	}
	if (dimcount != dimensions.size()) {
		LOGERROR("NetCDF variable %s/%s has %d dimensions, but data reference specifies %d dimensions", group.c_str(), variable_name.c_str(), dimcount, dimensions.size());
		return ref;
	}

	std::vector<std::string> mapped_indices(dimensions.size());
	for (size_t i = 0; i < dimcount; i++) {
		std::string dimname;
		if (!data_file.GetDimensionName(group, variable_name, i, dimname)) {
			return ref;
		} else {
			std::vector<std::string>::iterator it = std::find(dimensions.begin(), dimensions.end(), dimname);
			if (it == dimensions.end()) {
				LOGERROR("NetCDF variable %s/%s has a dimension %s which is not specified in the data reference dimensions", group.c_str(), variable_name.c_str(), dimname.c_str());
				return ref;
			} else {
				size_t offset = it - dimensions.begin();
				mapped_indices[i] = indices[offset];
			}
		}
	}

	ref = std::make_shared<DataReference>();
	if (!data_file.GetValue(group, variable_name, mapped_indices, &ref->value)) {
		ref.reset();
	}
	return ref;
}

}
