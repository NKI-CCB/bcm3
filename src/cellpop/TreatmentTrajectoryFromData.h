#pragma once

#include "TreatmentTrajectory.h"

class TreatmentTrajectoryFromData : public TreatmentTrajectory
{
public:
	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::NetCDFDataFile& data_file);
	virtual Real GetConcentration(Real time, Real cell_creation_time);

private:
	VectorReal timepoints;
	MatrixReal concentrations;
};
