#pragma once

class Experiment;

class TreatmentTrajectoryPulses : public TreatmentTrajectory
{
public:
	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::NetCDFDataFile& data_file);
	virtual Real GetConcentration(Real time, Real cell_creation_time);
	virtual Real FirstDiscontinuity(Real cell_creation_time);
	virtual Real NextDiscontinuity(Real time, Real cell_creation_time);

private:
	VectorReal timepoints;
};
