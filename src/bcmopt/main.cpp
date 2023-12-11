#include "Utils.h"
#include "LikelihoodFactory.h"
#include "Prior.h"
#include "ProgressIndicatorConsole.h"
#include "VariableSet.h"
#include "SampleHandlerStoreMaxAPosteriori.h"
#include "SamplerPT.h"
#include "NetCDFDataFile.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

int run(const po::variables_map& vm)
{
	std::shared_ptr<bcm3::VariableSet> varset;
	std::shared_ptr<bcm3::Prior> prior;
	std::shared_ptr<bcm3::Likelihood> likelihood;
	std::vector< std::shared_ptr<bcm3::Likelihood> > parallel_likelihoods;
	std::shared_ptr<bcm3::SamplerPT> sampler;
	std::shared_ptr<bcm3::SampleHandlerStoreMaxAPosteriori> output;

	std::string input_sample_file;
	size_t num_input_samples, start_sample_ix;
	std::string output_path;

	try {
		// Determine the number of threads we should use
		size_t numthreads = vm["sampling_threads"].as<size_t>();
		if (numthreads <= 0) {
			numthreads = std::thread::hardware_concurrency();
			printf("numthreads <= 0 specified, found %u hardware concurrency\n", (unsigned int)numthreads);
			if (numthreads <= 0) {
				numthreads = 1;
			}
		}
		size_t evaluation_threads = vm["evaluation_threads"].as<size_t>();

		varset = std::make_shared<bcm3::VariableSet>();
		if (!varset->LoadFromXML(vm["prior"].as<std::string>())) {
			return -2;
		}

		prior = bcm3::Prior::Create(vm["prior"].as<std::string>(), varset, numthreads);
		if (!prior) {
			return -3;
		}

		std::string likelihood_file = vm["likelihood"].as<std::string>();
		likelihood = bcm3::LikelihoodFactory::CreateLikelihood(likelihood_file, varset, vm, numthreads, evaluation_threads, true);
		if (!likelihood) {
			return -4;
		}
		if (!likelihood->IsReentrant()) {
			parallel_likelihoods.resize(numthreads);
			parallel_likelihoods[0] = likelihood;
			for (size_t i = 1; i < numthreads; i++) {
				parallel_likelihoods[i] = bcm3::LikelihoodFactory::CreateLikelihood(likelihood_file, varset, vm, numthreads, evaluation_threads, true);
			}
		}

		output_path = vm["output.folder"].as<std::string>();

		// Make sure output folder exists
		bcm3::fix_path(output_path);
		if (!boost::filesystem::is_directory(output_path)) {
			std::string create_path = output_path.substr(0, output_path.size() - 1);
			if (!boost::filesystem::create_directories(create_path)) {
				LOGERROR("Could not create output directory.");
				return false;
			}
		}

		// Open up a log file
		std::string log_filename = output_path + "log.txt";
		bcm3::logger->SetLogToFile(bcm3::Logger::Warning, log_filename.c_str());
		bcm3::logger->LogTime();

		// Create the sampler
		sampler = std::make_shared<bcm3::SamplerPT>(numthreads, 2ul * 1024 * 1024 * 1024);
		sampler->LoadSettings(vm);
		sampler->SetVariableSet(varset);
		sampler->SetPrior(prior);

		if (likelihood->IsReentrant()) {
			for (size_t i = 0; i < numthreads; i++) {
				sampler->SetLikelihood(likelihood, i);
			}
		} else {
			for (size_t i = 0; i < numthreads; i++) {
				sampler->SetLikelihood(parallel_likelihoods[i], i);
			}
		}

		output = std::make_shared<bcm3::SampleHandlerStoreMaxAPosteriori>();
		sampler->AddSampleHandler(output);
		sampler->SetOutputPath(output_path);

		// Everything has been set, initialize the sampler
		if (!sampler->Initialize()) {
			return -10;
		}

		input_sample_file = vm["input.sample_file"].as<std::string>();
		num_input_samples = vm["input.num_samples"].as<size_t>();
		start_sample_ix = vm["input.start_sample_ix"].as<size_t>();
	} catch (const std::exception & e) {
		LOGERROR("Program options error: %s", e.what());
		return -1;
	}

	//std::shared_ptr<bcm3::ProgressIndicator> pic = std::make_shared<bcm3::ProgressIndicatorConsole>();
	//sampler->SetProgressIndicator(pic);

	bcm3::NetCDFDataFile sample_file;
	if (!sample_file.Open(input_sample_file, false)) {
		LOGERROR("Unable to open input sample file \"%s\"", input_sample_file.c_str());
		return -11;
	}
	bool result = true;
	size_t total_samples, num_variables, num_temperatures;
	std::vector<std::string> variable_names;
	VectorReal temperatures;
	result &= sample_file.GetDimensionSize("samples", "sample_ix", &total_samples);
	result &= sample_file.GetDimensionSize("samples", "variable", &num_variables);
	result &= sample_file.GetDimensionSize("samples", "temperature", &num_temperatures);
	result &= sample_file.GetValues("samples", "variable", 0, num_variables, variable_names);
	result &= sample_file.GetValues("samples", "temperature", 0, num_temperatures, temperatures);

	std::vector<size_t> use_parameter_ix;
	std::vector<std::string> use_parameter_names;
	std::vector<unsigned int> use_variable_transforms;
	for (size_t i = 0; i < num_variables; i++) {
		bool is_sampled = false;
		for (size_t j = 0; j < varset->GetNumVariables(); j++) {
			if (varset->GetVariableName(j) == variable_names[i]) {
				is_sampled = true;
				break;
			}
		}
		if (!is_sampled) {
			use_parameter_ix.push_back(i);
			use_parameter_names.push_back(variable_names[i]);
			unsigned int transform;
			result &= sample_file.GetValue("samples", "variable_transform", i, &transform);
			use_variable_transforms.push_back(transform);
		}
	}

	if (!result) {
		return -12;
	}

	if (likelihood->IsReentrant()) {
		if (!likelihood->AddNonSampledParameters(use_parameter_names)) {
			return -7;
		}
		if (!likelihood->PostInitialize()) {
			return -6;
		}
	} else {
		for (size_t i = 0; i < parallel_likelihoods.size(); i++) {
			if (!parallel_likelihoods[i]->AddNonSampledParameters(use_parameter_names)) {
				return -7;
			}
			if (!parallel_likelihoods[i]->PostInitialize()) {
				return -6;
			}
		}
	}

	bcm3::logger->SetLogToConsole(bcm3::Logger::None);
	printf("Starting optimizations...\n");

	std::string outputfn = output_path + "MAP_estimates.tsv";
	FILE* file = fopen(outputfn.c_str(), "w");
	if (!file) {
		return -14;
	} else {
		fprintf(file, "temperature");
		for (size_t i = 0; i < num_input_samples; i++) {
			fprintf(file, "\t%zd", i);
		}
		fprintf(file, "\n");
		fclose(file);
	}

	std::string paramoutputfn = output_path + "MAP_estimates_paramvalues.tsv";
	file = fopen(paramoutputfn.c_str(), "w");
	if (!file) {
		return -14;
	} else {
		fprintf(file, "temperature_sample");
		for (size_t i = 0; i < use_parameter_ix.size(); i++) {
			fprintf(file, "\tfixed_%s", use_parameter_names[i].c_str());
		}
		for (size_t i = 0; i < varset->GetNumVariables(); i++) {
			fprintf(file, "\toptimized_%s", varset->GetVariableName(i).c_str());
		}
		fprintf(file, "\n");
		fclose(file);
	}

	for (size_t tempi = 0; tempi < num_temperatures; tempi++) {
		printf("Temperature %zd (%.6g)...\n", tempi, temperatures(tempi));
		printf("  Sample     ");
		for (size_t i = 0; i < num_input_samples; i++) {
			size_t sample_ix = start_sample_ix + i * (total_samples - start_sample_ix) / num_input_samples + ((total_samples - start_sample_ix) / num_input_samples - 1);
			printf("\b\b\b\b%4zd", i);
			fflush(stdout);

			VectorReal sample;
			if (!sample_file.GetValuesDim3("samples", "variable_values", sample_ix, tempi, 0, num_variables, sample)) {
				return -13;
			}
			VectorReal use_sample(use_parameter_ix.size());
			for (size_t j = 0; j < use_parameter_ix.size(); j++) {
				Real x = sample(use_parameter_ix[j]);
				if (use_variable_transforms[j] == bcm3::VariableSet::Transform_None) {
					use_sample(j) = x;
				} else if (use_variable_transforms[j] == bcm3::VariableSet::Transform_Log10) {
					use_sample(j) = bcm3::fastpow10(x);
				} else {
					LOGERROR("Not implemented");
					return -13;
				}
			}

			if (likelihood->IsReentrant()) {
				likelihood->SetNonSampledParameters(use_sample);
			} else {
				for (size_t j = 0; j < parallel_likelihoods.size(); j++) {
					parallel_likelihoods[j]->SetNonSampledParameters(use_sample);
				}
			}

			output->Reset();
			bool result = sampler->Run();
			Real MAP_lposterior;
			if (!result) {
				MAP_lposterior = std::numeric_limits<Real>::quiet_NaN();
			} else {
				MAP_lposterior = output->GetMAPlposterior();
			}

			file = fopen(outputfn.c_str(), "a");
			if (file) {
				if (i == 0) {
					fprintf(file, "%g", temperatures(tempi));
				}
				fprintf(file, "\t%g", MAP_lposterior);
				if (i == num_input_samples - 1) {
					fprintf(file, "\n");
				}
				fclose(file);
			} else {
				LOGERROR("Unable to open output file %s", outputfn.c_str());
			}

			file = fopen(paramoutputfn.c_str(), "a");
			if (file) {
				fprintf(file, "%g_%zd", temperatures(tempi), sample_ix);
				for (size_t pi = 0; pi < use_parameter_ix.size(); pi++) {
					fprintf(file, "\t%g", sample(use_parameter_ix[pi]));
				}
				for (size_t pi = 0; pi < varset->GetNumVariables(); pi++) {
					if (result) {
						fprintf(file, "\t%g", output->GetMAP()(pi));
					} else {
						fprintf(file, "\t%g", std::numeric_limits<Real>::quiet_NaN());
					}
				}
				fprintf(file, "\n");
				fclose(file);
			} else {
				LOGERROR("Unable to open output file %s", outputfn.c_str());
			}
		}
		printf("\n");
	}

	printf("Optimization finished.\n");
	bcm3::logger->SetLogToConsole(bcm3::Logger::Error);

	return 0;
}

int main(int argc, char* argv[])
{
	printf("BCM inference tool - version %u.%u.%u\n", BCM_VERSION_MAJOR, BCM_VERSION_MINOR, BCM_VERSION_FIX);

	bcm3::logger->SetLogToConsole(bcm3::Logger::Error);

	int result = -1;
	try {
		po::options_description generic_options("Generic options", 120);
		generic_options.add_options()
			("help", "produce help message")
			;

		po::options_description config_options("Configuration", 120);
		config_options.add_options()
			("sampling_threads,j",		po::value<size_t>()->default_value(0),						"number of sampling threads to use; set to 0 for using all available hardware concurrency")
			("evaluation_threads,k",	po::value<size_t>()->default_value(1),						"number of threads to use during each likelihood evaluation")
			("prior",					po::value<std::string>()->default_value("prior.xml"),		"File describing the priors for the parameters/initial conditions")
			("likelihood",				po::value<std::string>()->default_value("likelihood.xml"),	"File describing the likelihood functions")
			("output.folder",			po::value<std::string>()->default_value("output"),			"Folder in which output files will be generated, subfolder of the model directory.")
			("input.sample_file",		po::value<std::string>()->default_value("output.nc"),		"NetCDF file from which to load samples.")
			("input.num_samples",		po::value<size_t>()->default_value(50),						"Number of samples to select from the input sample file.")
			("input.start_sample_ix",	po::value<size_t>()->default_value(0),						"Start from this sample index.")
			;
#if 0
			("predict", "Make a prediction for the function described in the likelihood file, given a set of output samples.")
			("predict.input", po::value<std::string>()->default_value("output.nc"), "This should point to the output.nc file of a previous inference.")
			("predict.output", po::value<std::string>()->default_value("prediction.nc"), "Prediction output will be written to this file (in the output folder).")
			("tempered_prediction", po::value<std::string>(), "Option for generating predictions from tempered output samples.")
			("learning_rate,e", po::value<Real>()->default_value(1.0), "Learning rate between 0 and 1 - raises the likelihood function to this power.")
			("bootstrap_celllines", po::value<unsigned long>()->default_value(0), "Do inference with a bootstrap over cell lines")
			("safe_bayes_ix", po::value<size_t>()->default_value(std::numeric_limits<size_t>::max()), "SafeBayes - include cell lines up to and including index")
			("safe_bayes_cv_ix", po::value<std::string>()->default_value(std::string()), "Cross-validated SafeBayes - include cell lines excluding index")
			;
		//bcm::SamplerFactory::AddOptionsDescription(config_options);
#endif
		bcm3::SamplerPT::AddOptionsDescription(config_options);
		bcm3::LikelihoodFactory::AddOptionsDescription(config_options);
		po::options_description cmdline_only_options("Command-line only options", 120);
		cmdline_only_options.add_options()
			("config_file,c", po::value<std::string>()->default_value("config.txt"), "Configuration file")
			;

		// Add all options together
		po::options_description cmdline_options("Allowed options", 120);
		cmdline_options.add(generic_options).add(cmdline_only_options).add(config_options);

		// Parse command line
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);

		// Handle commands
		if (vm.count("help")) {
			std::cout << cmdline_options << "\n";
			result = 0;
		} else {
			// Parse config file, which may have been specified on the command line, otherwise load the default.
			std::string config_file = vm["config_file"].as<std::string>();
			if (!config_file.empty()) {
				po::store(po::parse_config_file<char>(config_file.c_str(), config_options), vm);
			}

#if 0
			if (vm.count("gen_variable_file")) {
				result = gen_variable_file(vm);
			} else if (vm.count("predict")) {
				result = predict(vm);
			} else {
				result = run(vm);
			}
#else
			result = run(vm);
#endif
		}
	} catch (po::error & e) {
		LOGERROR("Error processing program options: %s", e.what());
	}

	return result;
}
