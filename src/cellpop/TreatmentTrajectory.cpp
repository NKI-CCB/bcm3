#include "Utils.h"
#include "Experiment.h"
#include "TreatmentTrajectory.h"
#include "TreatmentTrajectoryFromData.h"
#include "TreatmentTrajectoryPulses.h"

std::unique_ptr<TreatmentTrajectory> TreatmentTrajectory::Create(const boost::property_tree::ptree& xml_node)
{
	std::string type = xml_node.get<std::string>("<xmlattr>.type");
	std::unique_ptr<TreatmentTrajectory> traj;
	if (type == "from_data") {
		traj = std::make_unique<TreatmentTrajectoryFromData>();
	} else if (type == "pulses") {
		traj = std::make_unique<TreatmentTrajectoryPulses>();
	} else {
		LOGERROR("Unknown trajectory type \"%s\"", type.c_str());
	}
	return traj;
}

bool TreatmentTrajectory::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::NetCDFDataFile& data_file)
{
	return true;
}

Real TreatmentTrajectory::FirstDiscontinuity(Real cell_creation_time)
{
	return std::numeric_limits<OdeReal>::quiet_NaN();
}

Real TreatmentTrajectory::NextDiscontinuity(Real time, Real cell_creation_time)
{
	return std::numeric_limits<OdeReal>::quiet_NaN();
}
