#pragma once

inline Real GetDrugMolecularWeight(const std::string& drug)
{
	if (drug == "lapatinib") {
		return 581.06;
	} else if (drug == "dacomitinib") {
		return 469.95;
	} else if (drug == "afatinib") {
		return 485.94;
	} else if (drug == "trametinib") {
		return 615.404;
	} else if (drug == "mirdametinib") {
		return 482.19;
	} else if (drug == "selumetinib") {
		return 457.68;
	} else {
		LOGERROR("Unknown drug \"%s\"", drug.c_str());
		return std::numeric_limits<Real>::quiet_NaN();
	}
}
