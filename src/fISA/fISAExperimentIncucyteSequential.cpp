#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "fISAExperimentIncucyteSequential.h"
#include "fISAExperimentSingleCondition.h"
#include "SignalingNetwork.h"
#include "VectorUtils.h"

#include <boost/property_tree/xml_parser.hpp>
#include "CSVParser.h"

#define USE_VINE_COPULAS 0
#define USE_MVT 0

fISAExperimentIncucyteSequential::fISAExperimentIncucyteSequential()
	: relative_experiment(NULL)
	, drug_model_ix(std::numeric_limits<size_t>::max())
{
}

fISAExperimentIncucyteSequential::~fISAExperimentIncucyteSequential()
{
}

bool fISAExperimentIncucyteSequential::StartEvaluateLogProbability(const VectorReal& values, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	for (size_t vi = 0; vi < varset->GetNumVariables(); vi++) {
		transformed_values[vi] = varset->TransformVariable(vi, values[vi]);
		if (transformed_values[vi] == std::numeric_limits<Real>::infinity()) {
			LOGWARNING("Overflow of parameter %u after transformation of value %g; truncating to highest floating point value", vi, values[vi]);
			transformed_values[vi] = std::numeric_limits<Real>::max();
		} else if (transformed_values[vi] == -std::numeric_limits<Real>::infinity()) {
			LOGWARNING("Negative overflow of parameter %u after transformation of value %g; truncating to lowest floating point value", vi, values[vi]);
			transformed_values[vi] = std::numeric_limits<Real>::lowest();
		}
	}

	for (ptrdiff_t ci = 0; ci < cell_lines.size(); ci++) {
		bcm3::TTask task = boost::bind(&fISAExperimentIncucyteSequential::EvaluateCellLine, this, boost::placeholders::_1, boost::placeholders::_2);
		evaluation_tasks[ci] = task_manager->AddTask(task, (void*)ci);
	}

	return true;
}

bool fISAExperimentIncucyteSequential::FinishEvaluateLogProbability(const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	logp = 0.0;
	for (ptrdiff_t ci = 0; ci < cell_lines.size(); ci++) {
		bool result = task_manager->WaitTask(evaluation_tasks[ci]);
		if (result) {
			logp += cell_line_logp(ci);
		} else {
			LOGERROR("Evaluation %u failed", ci);
			return false;
		}
	}

	return true;
}

void fISAExperimentIncucyteSequential::GetObservedData(size_t data_ix, MatrixReal& out_values) const
{
	out_values.resize(cell_lines.size(), 1);
	for (ptrdiff_t ci = 0; ci < cell_lines.size(); ci++) {
		size_t dci = data_ix / 2;
		if (data_ix % 2 == 0) {
			out_values(ci, 0) = data.gmms[ci * drug_concentrations.size() + dci].mup[0];
		} else {
			out_values(ci, 0) = data.gmms[ci * drug_concentrations.size() + dci].mua[0];
		}
	}
}

void fISAExperimentIncucyteSequential::GetModeledData(size_t threadix, size_t data_ix, VectorReal& out_values) const
{
	out_values.resize(cell_lines.size());
	for (ptrdiff_t ci = 0; ci < cell_lines.size(); ci++) {
		size_t dci = data_ix / 2;
		if (data_ix % 2 == 0) {
			out_values(ci) = modeled_data_values[ci][dci](0);
		} else {
			out_values(ci) = modeled_data_values[ci][dci](1);
		}
	}
}

void fISAExperimentIncucyteSequential::GetModeledActivities(size_t threadix, MatrixReal& out_values) const
{
	out_values.resize(cell_lines.size(), network->GetMoleculeCount());
	for (ptrdiff_t ci = 0; ci < cell_lines.size(); ci++) {
		out_values.row(ci) = activities[ci][0];
	}
}

bool fISAExperimentIncucyteSequential::LoadTypeSpecificNodes(const boost::property_tree::ptree& xml_node, const bcm3::NetCDFDataFile& datafile)
{
	bool result = true;

	try {
		const boost::property_tree::ptree& drug_range = xml_node.get_child("drug_range");

		drug_species_name = drug_range.get<std::string>("<xmlattr>.species_name");
		drug_model_ix = network->GetSignalingMoleculeIxByName(drug_species_name);
		if (drug_model_ix == std::numeric_limits<size_t>::max()) {
			LOGERROR("Could not find signaling node \"%s\" for drug range data", drug_species_name.c_str());
			return false;
		}

		std::string concentrations = drug_range.get<std::string>("<xmlattr>.concentrations", "");
		if (concentrations.empty()) {
			std::string data_name = drug_range.get<std::string>("<xmlattr>.concentrations_data_name");

			size_t concentration_count = 0;
			result &= datafile.GetDimensionSize(Name, data_name, &concentration_count);
			result &= datafile.GetValues(Name, data_name, 0, concentration_count, drug_concentrations);
		} else {
			if (!bcm3::ParseVectorFromString(concentrations, drug_concentrations)) {
				LOGERROR("Unable to parse concentration string \"%s\"", concentrations.c_str());
				return false;
			}
		}
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		return false;
	}

	return result;
}

bool fISAExperimentIncucyteSequential::InitializeParallelData(size_t evaluation_threads)
{
	parallel_data.resize(evaluation_threads);

	cell_line_logp.setZero(cell_lines.size());
	evaluation_tasks.resize(cell_lines.size(), 0);
	activities.resize(cell_lines.size());
	modeled_data_values.resize(cell_lines.size());
	for (ptrdiff_t ci = 0; ci < cell_lines.size(); ci++) {
		activities[ci].resize(drug_concentrations.size());
		modeled_data_values[ci].resize(drug_concentrations.size());
		for (ptrdiff_t dci = 0; dci < drug_concentrations.size(); dci++) {
			activities[ci][dci] = VectorReal::Constant(network->GetMoleculeCount(), std::numeric_limits<Real>::quiet_NaN());
			modeled_data_values[ci][dci] = VectorReal::Constant(2, std::numeric_limits<Real>::quiet_NaN());
		}
	}

	return true;
}

bool fISAExperimentIncucyteSequential::ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile)
{
	bool result = true;
	
	std::string data_file_base = xml_node.get<std::string>("<xmlattr>.data_file_base");

	size_t cell_line_count;
	result &= datafile.GetDimensionSize(Name, "cell_lines", &cell_line_count);

	data.model_ix = network->GetSignalingMoleculeIxByName("proliferation");
	data.model_ix2 = network->GetSignalingMoleculeIxByName("apoptosis");

#if USE_VINE_COPULAS
	data.vinecops.resize(cell_line_count * drug_concentrations.size(), bcm::VineCopula(2));
	boost::property_tree::ptree pt;
	try {
		boost::property_tree::read_xml(data_file_base, pt);
		boost::property_tree::ptree prior_node = pt.get_child("response_estimate");

		BOOST_FOREACH(const boost::property_tree::ptree::value_type& var, prior_node.get_child("")) {
			if (var.first == "prior") {
				std::string cell_line = var.second.get<std::string>("<xmlattr>.cell_line");
				size_t ci = var.second.get<size_t>("<xmlattr>.concentration") - 1;
				size_t i;
				result &= datafile->GetDimensionIx(Name, "cell_lines", cell_line, &i);
				result &= data.vinecops[i * 9 + ci].LoadFromXML(var.second);
			}
		}

	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing variable file: %s", e.what());
		return NULL;
	}
#elif USE_MVT
	data.tdists.resize(cell_line_count * drug_concentrations.size());
	bcm3::CSVParser parser;
	parser.Parse(data_file_base, "\t");

	for (size_t i = 0; i < cell_line_count; i++) {
		std::string cell_line;
		result &= datafile->GetValue(Name, "cell_lines", i, cell_line);

		for (int ci = 0; ci < drug_concentrations.size(); ci++) {
			data.tdists[i * 9 + ci].mup = parser.GetEntry(i * 9 + ci, 0);
			data.tdists[i * 9 + ci].mua = parser.GetEntry(i * 9 + ci, 1);
			data.tdists[i * 9 + ci].cov = MatrixReal::Zero(2, 2);
			data.tdists[i * 9 + ci].cov(0, 0) = parser.GetEntry(i * 9 + ci, 2);
			data.tdists[i * 9 + ci].cov(0, 1) = parser.GetEntry(i * 9 + ci, 3);
			data.tdists[i * 9 + ci].cov(1, 0) = parser.GetEntry(i * 9 + ci, 3);
			data.tdists[i * 9 + ci].cov(1, 1) = parser.GetEntry(i * 9 + ci, 4);
			data.tdists[i * 9 + ci].invcov = data.tdists[i * 9 + ci].cov.inverse();
		}
	}
#else
	data.gmms.resize(cell_line_count * drug_concentrations.size());
	bcm3::CSVParser parser;
	parser.Parse(data_file_base, "\t");

	for (ptrdiff_t i = 0; i < cell_line_count; i++) {
		std::string cell_line;
		result &= datafile.GetValue(Name, "cell_lines", i, cell_line);

		for (ptrdiff_t ci = 0; ci < drug_concentrations.size(); ci++) {
			for (int ki = 0; ki < 3; ki++) {
				DataPart::gmm& gm = data.gmms[i * drug_concentrations.size() + ci];

				gm.mup[ki] = parser.GetEntry(i * 9 + ci, ki * 5 + 0);
				gm.mua[ki] = parser.GetEntry(i * 9 + ci, ki * 5 + 1);
				gm.cov[ki] = MatrixReal::Zero(2, 2);
				gm.cov[ki](0, 0) = parser.GetEntry(i * 9 + ci, ki * 5 + 2);
				gm.cov[ki](0, 1) = parser.GetEntry(i * 9 + ci, ki * 5 + 3);
				gm.cov[ki](1, 0) = parser.GetEntry(i * 9 + ci, ki * 5 + 3);
				gm.cov[ki](1, 1) = parser.GetEntry(i * 9 + ci, ki * 5 + 4);
				gm.invcov[ki] = gm.cov[ki].inverse();
				gm.weights[ki] = parser.GetEntry(i * 9 + ci, 15 + ki);
				gm.logncweight[ki] = log(gm.weights[ki]) - log(2 * M_PI * sqrt(gm.cov[ki].determinant()));
			}
		}
	}
#endif

	// Test for relative experiment.
	std::string type = xml_node.get<std::string>("<xmlattr>.type", "");
	if (type == "relative") {
		std::string relative_to_experiment = xml_node.get<std::string>("<xmlattr>.relative_reference");

		// Make sure the other experiment is defined before this one
		bool found = false;
		for (std::vector< std::unique_ptr<fISAExperiment> >::const_iterator ei = other_experiments.begin(); ei != other_experiments.end(); ++ei) {
			if (ei->get()->GetName() == relative_to_experiment) {
				relative_experiment = dynamic_cast<const fISAExperimentSingleCondition*>(ei->get());
				if (relative_experiment == NULL) {
					LOGERROR("Relative experiment \"%s\" is not a single condition experiment", relative_to_experiment.c_str());
					return false;
				} else {
					found = true;
					break;
				}
			}
		}
		if (!found) {
			LOGERROR("A data point for experiment \"%s\" is specified as relative to experiment \"%s\", but that experiment has not been defined or is not of type single condition. Make sure relative experiments are defined in the right order.", Name.c_str(), relative_to_experiment.c_str());
			return false;
		}
	}

	return result;
}

bool fISAExperimentIncucyteSequential::EvaluateCellLine(void* cell_line_ix_as_ptr, size_t eval_thread_ix)
{
	ParallelDataBase& pd = parallel_data[eval_thread_ix];

	size_t cell_line_ix = (size_t)cell_line_ix_as_ptr;

	Real logp = 0.0;
	for (ptrdiff_t dci = 0; dci < drug_concentrations.size(); dci++) {
		VectorReal& activities = this->activities[cell_line_ix][dci];

		// Calculate model activities
		PrepareActivitiesCalculation(activities, pd.expression, transformed_values.data(), cell_line_ix);
		activities(drug_model_ix) = drug_concentrations(dci);
		if (network->Calculate(eval_thread_ix, activities, pd.expression, transformed_values.data())) {
			Real proliferation = activities(data.model_ix);
			Real apoptosis = activities(data.model_ix2);

			modeled_data_values[cell_line_ix][dci](0) = proliferation;
			modeled_data_values[cell_line_ix][dci](1) = apoptosis;

			if (relative_experiment) {
				const Real relative_reference_proliferation = relative_experiment->stored_activities[cell_line_ix](data.model_ix);
				proliferation = proliferation - relative_reference_proliferation;
			}

			Real thisp;
			Real x[2] = { proliferation, apoptosis };
#if USE_VINE_COPULAS
			try {
				if (!data.vinecops[ci * 9 + dci].EvaluateLogPDF(x, thisp)) {
					logp = -std::numeric_limits<Real>::infinity();
					break;
				} else {
					logp += thisp;
				}
			} catch (std::exception&) {
				logp = -std::numeric_limits<Real>::infinity();
				break;
			}
#elif USE_MVT
			const DataPart::tdist& td = data.tdists[ci * 9 + dci];
			const Real lognc = log(sqrt(td.invcov.determinant()) / (2 * M_PI));
			const Real tx = proliferation - td.mup;
			const Real ta = apoptosis - td.mua;
			thisp = lognc -
				log1p((1.0 / 3.0) * td.invcov(0, 0) * tx * tx +
				(1.0 / 3.0) * td.invcov(1, 1) * ta * ta +
					(1.0 / 3.0) * td.invcov(0, 1) * tx * ta +
					(1.0 / 3.0) * td.invcov(1, 0) * tx * ta) *
				2.5;
			logp += thisp;
#else
			const DataPart::gmm& gm = data.gmms[cell_line_ix * drug_concentrations.size() + dci];
			if (gm.mua[1] == gm.mua[1] && gm.mup[1] == gm.mup[1]) {
				thisp = -std::numeric_limits<Real>::infinity();
				for (int ki = 0; ki < 3; ki++) {
					const Real tx = proliferation - gm.mup[ki];
					const Real ta = apoptosis - gm.mua[ki];
					//const Real kp = gm.logncweight[ki] - 0.5 * (gm.invcov[ki](0, 0) * tx * tx +
					//											gm.invcov[ki](1, 1) * ta * ta +
					//											gm.invcov[ki](1, 0) * ta * tx +
					//											gm.invcov[ki](0, 1) * tx * ta );
					if (gm.weights[ki] > 0) {
						const Real kp = gm.logncweight[ki] - 2.5 * log1p((1.0 / 3.0) * gm.invcov[ki](0, 0) * tx * tx +
							(1.0 / 3.0) * gm.invcov[ki](1, 1) * ta * ta +
							(1.0 / 3.0) * gm.invcov[ki](1, 0) * ta * tx +
							(1.0 / 3.0) * gm.invcov[ki](0, 1) * tx * ta);
						thisp = bcm3::logsum(thisp, kp);
					}
				}
				logp += thisp;
			}
#endif
		} else {
			logp = -std::numeric_limits<Real>::infinity();
			break;
		}
	}

	cell_line_logp(cell_line_ix) = logp;

	return true;
}
