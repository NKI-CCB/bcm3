#pragma once

#include "NetCDFDataFile.h"

class Experiment;

class TreatmentTrajectory
{
	// Treatment trajectories are defined in global time, but the cell specific time is provided
	// to ensure exact floating point times for the discontinuities
public:
	static std::unique_ptr<TreatmentTrajectory> Create(const boost::property_tree::ptree& xml_node);

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::NetCDFDataFile& data_file);
	virtual Real GetConcentration(Real time, Real cell_creation_time) = 0;
	virtual Real FirstDiscontinuity(Real cell_creation_time);
	virtual Real NextDiscontinuity(Real time, Real cell_creation_time);
};
