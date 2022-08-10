#pragma once

namespace bcm3 {

class CSVParser
{
public:
	bool Parse(const std::string& filename, const char* sep = ",", bool rownames = true);

	inline size_t GetNumColumns() const { return ColumnNames.size(); }
	inline size_t GetNumRows() const { return Rows.size(); }
	inline const std::string& GetColumnName(size_t col) const { return ColumnNames[col]; }
	inline const std::string& GetRowName(size_t row) const { return RowNames[row]; }
	inline Real GetEntry(size_t row, size_t col) const { return Rows[row][col]; }

private:
	std::vector<std::string> ColumnNames;
	std::vector<std::string> RowNames;
	typedef std::vector<Real> Row;
	std::vector<Row> Rows;
};

}
