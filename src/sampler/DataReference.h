#pragma once

#include "NetCDFDataFile.h"

namespace bcm3 {

class DataReference
{
public:
	DataReference();
	~DataReference();

	static std::shared_ptr<DataReference> Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::string dim1, std::string ix1);
	static std::shared_ptr<DataReference> Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::string dim1, std::string ix1, std::string dim2, std::string ix2);
	static std::shared_ptr<DataReference> Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::string dim1, std::string ix1, std::string dim2, std::string ix2, std::string dim3, std::string ix3);
	static std::shared_ptr<DataReference> Create(const NetCDFDataFile& data_file, std::string group, std::string variable_name, std::vector<std::string> dimensions, std::vector<std::string> indices);

	inline Real GetValue() const { return value; }

private:
	Real value;
};

}
