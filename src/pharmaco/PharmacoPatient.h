#pragma once

#include "NetCDFDataFile.h"

struct Patient
{
	Patient();

	bool Load(const bcm3::NetCDFDataFile& data, const std::string& trial, const std::string& drug);

	std::string patient_id;

	VectorReal treatment_timepoints;
	VectorReal treatment_doses;
	VectorReal observation_timepoints;
	VectorReal observed_concentrations;
	VectorReal simulated_concentrations;
};
