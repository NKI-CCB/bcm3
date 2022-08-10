#include "Utils.h"
#include "CSVParser.h"

#include <fstream>
#include <boost/algorithm/string/trim.hpp>

namespace bcm3 {

bool CSVParser::Parse(const std::string& filename, const char* sep, bool rownames)
{
	std::ifstream ifs;
	ifs.open(filename.c_str(), std::ifstream::in);
	if (!ifs.is_open()) {
		LOGERROR("Could not open file \"%s\" for reading", filename.c_str());
		return false;
	}

	std::string line;
	if (!getline(ifs, line)) {
		LOGERROR("Could not read from file \"%s\"", filename.c_str());
		ifs.close();
		return false;
	}
	if (*line.rbegin() == '\r') {
		line.pop_back();
	}
	
	bcm3::tokenize(line, ColumnNames, sep);
	if (ColumnNames.empty()) {
		LOGERROR("CSV file \"%s\" has no header with column names", filename.c_str());
		ifs.close();
		return false;
	}
	if (rownames) {
		ColumnNames.erase(ColumnNames.begin());
	}

	size_t data_line = 0;
	while(getline(ifs, line)) {
		if (*line.rbegin() == '\r') {
			line.pop_back();
		}
		std::vector<std::string> rowstr;
		rowstr.reserve(ColumnNames.size());
		bcm3::tokenize(line, rowstr, sep);

		if ( (rownames && rowstr.size() != ColumnNames.size() + 1) ||
			(!rownames && rowstr.size() != ColumnNames.size()    )) {
			LOGERROR("Inconsistent CSV data file \"%s\": found %d instead of %d data points at data line %d", filename.c_str(), rowstr.size(), ColumnNames.size(), data_line);
			return false;
		}

		if (rownames) {
			RowNames.push_back(rowstr[0]);
		} else {
			char tmp[16];
			snprintf(tmp, 16, "%llu", (uint64)data_line);
			RowNames.push_back(tmp);
		}

		Row row;
		row.resize(rownames ? rowstr.size() - 1 : rowstr.size());
		for (size_t i = 0; i < ColumnNames.size(); i++) {
			const char* value = rownames ? rowstr[i+1].c_str() : rowstr[i].c_str();
			if (strcasecmp(value, "NA") == 0 || strcasecmp(value, "nan") == 0) {
				row[i] = std::numeric_limits<Real>::quiet_NaN();
			} else {
				row[i] = atof(value);
			}
		}

		Rows.push_back(row);
		data_line++;
	}

	ifs.close();
	return true;
}

}
