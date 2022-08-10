#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "fISAExperimentIncucyteSequential.h"
#include "fISAExperimentSingleCondition.h"
#include "SignalingNetwork.h"
#include "VectorUtils.h"

#include <boost/property_tree/xml_parser.hpp>
#include "CSVParser.h"

#if TODO

#define USE_VINE_COPULAS 0
#define USE_MVT 0

fISAExperimentIncucyteSequential::fISAExperimentIncucyteSequential(size_t numthreads)
	: fISAExperiment(numthreads)
	, relative_experiment(NULL)
	, data_weight(1.0)
{
}

fISAExperimentIncucyteSequential::~fISAExperimentIncucyteSequential()
{
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

bool fISAExperimentIncucyteSequential::InitializeParallelData()
{
	parallel_data.resize(numthreads);
	for (size_t threadix = 0; threadix < numthreads; threadix++) {
		parallel_data[threadix].transformed_values.resize(varset->GetNumVariables(), std::numeric_limits<Real>::quiet_NaN());
		parallel_data[threadix].activities.resize(cell_lines.size());
		parallel_data[threadix].modeled_data_values.resize(cell_lines.size());
		for (size_t ci = 0; ci < cell_lines.size(); ci++) {
			parallel_data[threadix].activities[ci].resize(drug_concentrations.size(), VectorReal::Constant(network->GetMoleculeCount(), std::numeric_limits<Real>::quiet_NaN()));
			parallel_data[threadix].modeled_data_values[ci].resize(drug_concentrations.size(), VectorReal::Constant(2, std::numeric_limits<Real>::quiet_NaN()));
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

	for (size_t i = 0; i < cell_line_count; i++) {
		std::string cell_line;
		result &= datafile.GetValue(Name, "cell_lines", i, cell_line);

		for (int ci = 0; ci < drug_concentrations.size(); ci++) {
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
			if (ei->get()->GetName() == relative_to_experiment && (relative_experiment = dynamic_cast<fISAExperimentSingleCondition*>(ei->get())) != NULL) {
				found = true;
				break;
			}
		}
		if (!found) {
			LOGERROR("A data point for experiment \"%s\" is specified as relative to experiment \"%s\", but that experiment has not been defined or is not of type single condition. Make sure relative experiments are defined in the right order.", Name.c_str(), relative_to_experiment.c_str());
			return false;
		}
	}

	data_weight = xml_node.get<Real>("<xmlattr>.weight", 1.0);

	return result;
}

bool fISAExperimentIncucyteSequential::EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments)
{
	ParallelData& pd = parallel_data[threadix];

	for (size_t vi = 0; vi < varset->GetNumVariables(); vi++) {
		pd.transformed_values[vi] = varset->TransformVariable(vi, values[vi]);
	}

	logp = 0.0;
	for (size_t ci = 0; ci < cell_lines.size(); ci++) {
		for (int dci = 0; dci < drug_concentrations.size(); dci++) {
			VectorReal& activities = pd.activities[ci][dci];

			// Calculate model activities
			PrepareActivitiesCalculation(activities, pd.expression, pd.transformed_values.data(), ci);
			activities(drug_model_ix) = drug_concentrations(dci);
			if (network->Calculate(threadix, activities, pd.expression, pd.transformed_values.data())) {
				Real proliferation = activities(data.model_ix);
				Real apoptosis = activities(data.model_ix2);

				pd.modeled_data_values[ci][dci](0) = proliferation;
				pd.modeled_data_values[ci][dci](1) = apoptosis;

#if 0
				if (relative_experiment) {
					const Real relative_reference_proliferation = relative_experiment->parallel_data[threadix].activities[ci](data.model_ix);
					proliferation = proliferation - relative_reference_proliferation;
				}
#endif

				Real thisp;
				Real x[2] = { proliferation, apoptosis };
#if USE_VINE_COPULAS
				try {
					if (!data.vinecops[ci * 9 + dci].EvaluateLogPDF(x, thisp)) {
						logp = -std::numeric_limits<Real>::infinity();
						break;
					} else {
						logp += data_weight * thisp;
					}
				} catch (std::exception&) {
					logp = -std::numeric_limits<Real>::infinity();
					break;
				}
#elif USE_MVT
				const DataPart::tdist& td =	data.tdists[ci * 9 + dci];
				const Real lognc = log(sqrt(td.invcov.determinant()) / (2 * M_PI));
				const Real tx = proliferation - td.mup;
				const Real ta = apoptosis - td.mua;
				thisp = lognc -
						log1p(	(1.0/3.0) * td.invcov(0, 0) * tx * tx + 
								(1.0/3.0) * td.invcov(1, 1) * ta * ta +
								(1.0/3.0) * td.invcov(0, 1) * tx * ta +
								(1.0/3.0) * td.invcov(1, 0) * tx * ta ) *
						2.5;
				logp += data_weight * thisp;
#else
				const DataPart::gmm& gm = data.gmms[ci * drug_concentrations.size() + dci];
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
					logp += data_weight * (1.0 / drug_concentrations.size()) * thisp;
				}
#endif
			} else {
				logp = -std::numeric_limits<Real>::infinity();
				break;
			}
		}

		if (logp == -std::numeric_limits<Real>::infinity()) {
			break;
		}
	}

	return true;
}

#endif
