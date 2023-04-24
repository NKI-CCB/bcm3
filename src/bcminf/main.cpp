#include "Utils.h"
#include "LikelihoodFactory.h"
#include "Prior.h"
#include "ProgressIndicatorConsole.h"
#include "VariableSet.h"
#include "SampleHandlerNetCDF.h"
#include "SamplerFactory.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

int run(const po::variables_map& vm)
{
	std::shared_ptr<bcm3::VariableSet> varset;
	std::shared_ptr<bcm3::Prior> prior;
	std::shared_ptr<bcm3::Likelihood> likelihood;
	std::vector< std::shared_ptr<bcm3::Likelihood> > parallel_likelihoods;
	std::shared_ptr<bcm3::Sampler> sampler;
	std::shared_ptr<bcm3::SampleHandlerNetCDF> output;
	std::string output_path;
	double progress_update_time = 0.5;

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
		likelihood = bcm3::LikelihoodFactory::CreateLikelihood(likelihood_file, varset, vm, numthreads, evaluation_threads);
		if (!likelihood) {
			return -4;
		}
		if (!likelihood->IsReentrant()) {
			parallel_likelihoods.resize(numthreads);
			parallel_likelihoods[0] = likelihood;
			for (size_t i = 1; i < numthreads; i++) {
				parallel_likelihoods[i] = bcm3::LikelihoodFactory::CreateLikelihood(likelihood_file, varset, vm, numthreads, evaluation_threads);
			}
		}

		progress_update_time = vm["progress_update_time"].as<double>();
		output_path = vm["output.folder"].as<std::string>();

		// Make sure output folder exists
		bcm3::fix_path(output_path);
		if (!boost::filesystem::is_directory(output_path)) {
			std::string create_path = output_path.substr(0, output_path.size() - 1);
			if (!boost::filesystem::create_directories(create_path)) {
				LOGERROR("Could not create output directory.");
				return -5;
			}
		}

		// Open up a log file
		std::string log_filename = output_path + "log.txt";
		bcm3::logger->SetLogToFile(bcm3::Logger::Info, log_filename.c_str());
		bcm3::logger->LogTime();

		// Create the sampler
		sampler = bcm3::SamplerFactory::Create(vm, numthreads, 2ul * 1024 * 1024 * 1024);
		sampler->SetVariableSet(varset);
		sampler->SetPrior(prior);

		Real learning_rate = vm["learning_rate"].as<Real>();
		if (likelihood->IsReentrant()) {
			for (size_t i = 0; i < numthreads; i++) {
				likelihood->SetLearningRate(learning_rate);
				sampler->SetLikelihood(likelihood, i);
			}
		} else {
			for (size_t i = 0; i < numthreads; i++) {
				parallel_likelihoods[i]->SetLearningRate(learning_rate);
				sampler->SetLikelihood(parallel_likelihoods[i], i);
			}
		}

		output = std::make_shared<bcm3::SampleHandlerNetCDF>();
		output->SetFile(output_path + "output.nc");
		output->Initialize(sampler->GetNumOutputSamples(), varset.get(), sampler->GetOutputTemperatures());
		sampler->AddSampleHandler(output);
		sampler->SetOutputPath(output_path);

		// Everything has been set, initialize the sampler
		if (!sampler->Initialize()) {
			return -10;
		}
	} catch (const std::exception & e) {
		LOGERROR("Program options error: %s", e.what());
		return -1;
	}

	if (likelihood->IsReentrant()) {
		if (!likelihood->PostInitialize()) {
			return -6;
		}
	} else {
		for (size_t i = 0; i < parallel_likelihoods.size(); i++) {
			if (!parallel_likelihoods[i]->PostInitialize()) {
				return -6;
			}
		}
	}

	std::shared_ptr<bcm3::ProgressIndicatorConsole> pic = std::make_shared<bcm3::ProgressIndicatorConsole>();
	pic->SetUpdateTime(progress_update_time);
	sampler->SetProgressIndicator(pic);

	printf("Starting sampler...\n");
	bool result = sampler->Run();
	printf("Sampling finished.\n");

	output->Close();

	likelihood->OutputEvaluationStatistics(output_path);

	if (result) {
		return 0;
	} else {
		return -20;
	}
}

int predict(const po::variables_map& vm)
{
	std::shared_ptr<bcm3::VariableSet> varset;
	std::shared_ptr<bcm3::Prior> prior;
	std::shared_ptr<bcm3::Likelihood> likelihood;
	std::vector< std::shared_ptr<bcm3::Likelihood> > parallel_likelihoods;

	std::string predict_input_fn;
	std::string predict_output_fn;

	size_t specific_temperature;
	size_t skip_n;

	try {
		size_t evaluation_threads = vm["evaluation_threads"].as<size_t>();
		varset = std::make_shared<bcm3::VariableSet>();
		if (!varset->LoadFromXML(vm["prior"].as<std::string>())) {
			return -2;
		}

		prior = bcm3::Prior::Create(vm["prior"].as<std::string>(), varset, 1);
		if (!prior) {
			return -3;
		}

		std::string likelihood_file = vm["likelihood"].as<std::string>();
		likelihood = bcm3::LikelihoodFactory::CreateLikelihood(likelihood_file, varset, vm, 1, evaluation_threads);
		if (!likelihood) {
			return -4;
		}

		specific_temperature = vm["predict.specific_temperature"].as<size_t>();
		skip_n = vm["predict.skip_n"].as<size_t>();

		std::string output_path = vm["output.folder"].as<std::string>();
		predict_input_fn = vm["predict.input"].as<std::string>();
		predict_output_fn = vm["predict.output"].as<std::string>();

		// Make sure output folder exists
		bcm3::fix_path(output_path);
		if (!boost::filesystem::is_directory(output_path)) {
			std::string create_path = output_path.substr(0, output_path.size() - 1);
			if (!boost::filesystem::create_directories(create_path)) {
				LOGERROR("Could not create output directory.");
				return -5;
			}
		}

		// Open up a log file
		std::string log_filename = output_path + "log.txt";
		bcm3::logger->SetLogToFile(bcm3::Logger::Info, log_filename.c_str());
		bcm3::logger->LogTime();

		predict_output_fn = output_path + predict_output_fn;
	} catch (const std::exception & e) {
		LOGERROR("Program options error: %s", e.what());
		return -1;
	}

	if (!likelihood->PostInitialize()) {
		return -6;
	}

	bcm3::NetCDFDataFile sample_file;
	if (!sample_file.Open(predict_input_fn, false)) {
		LOGERROR("Unable to open input sample file \"%s\"", predict_input_fn.c_str());
		return -11;
	}
	bool result = true;
	size_t total_samples, num_variables, num_temperatures;
	std::vector<std::string> variable_names;
	VectorReal temperatures;
	std::vector<unsigned int> sample_ixs;
	result &= sample_file.GetDimensionSize("samples", "sample_ix", &total_samples);
	result &= sample_file.GetDimensionSize("samples", "variable", &num_variables);
	result &= sample_file.GetDimensionSize("samples", "temperature", &num_temperatures);
	result &= sample_file.GetValues("samples", "sample_ix", 0, total_samples, sample_ixs);
	result &= sample_file.GetValues("samples", "variable", 0, num_variables, variable_names);
	result &= sample_file.GetValues("samples", "temperature", 0, num_temperatures, temperatures);
	if (!result) {
		return -12;
	}

	bcm3::NetCDFDataFile output_file;
	if (!output_file.Create(predict_output_fn)) {
		LOGERROR("Unable to open output file \"%s\"", predict_output_fn.c_str());
		return -13;
	}
	result &= output_file.CreateGroup("predictions");
	result &= output_file.CreateDimension("predictions", "sample_ix", sample_ixs);
	result &= output_file.CreateDimension("predictions", "temperature", temperatures);
	result &= output_file.CreateVariable("predictions", "log_likelihood", "sample_ix", "temperature");
	if (!result) {
		return -14;
	}

	bcm3::Timer timer;
	timer.Start();

	unsigned int num_likelihood_evaluations = 0;
	for (int j = 0; j < temperatures.size(); j++) {
		if (specific_temperature == std::numeric_limits<size_t>::max() || (int)specific_temperature == j) {
			printf("Temperature %d (%.6g)...\n", j, temperatures(j));
			printf("  Sample     ");
			for (size_t i = total_samples / 2; i < total_samples; i += skip_n + 1) {
				printf("\b\b\b\b%4zd", i);
				fflush(stdout);

				VectorReal sample;
				if (!sample_file.GetValuesDim3("samples", "variable_values", i, j, 0, num_variables, sample)) {
					return -12;
				}

				Real llh = 0.0;
				result = likelihood->EvaluateLogProbability(0, sample, llh);
				if (result) {
				} else {
					llh = std::numeric_limits<Real>::quiet_NaN();
				}
				num_likelihood_evaluations++;

				output_file.PutValue("predictions", "log_likelihood", i, j, llh);
			}
			printf("\n");
		}
	}

	double prediction_time = timer.GetElapsedSeconds();
	LOG("  Integration time: %.3g seconds", prediction_time);
	LOG("  Number of likelihood evaluations: %u", num_likelihood_evaluations);
	LOG("  Evaluations per second: %.0f", num_likelihood_evaluations / prediction_time);

	sample_file.Close();
	output_file.Close();

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
			("sampling_threads,j",				po::value<size_t>()->default_value(0),									"number of sampling threads to use; set to 0 for using all available hardware concurrency")
			("evaluation_threads,k",			po::value<size_t>()->default_value(1),									"number of threads to use during each likelihood evaluation")
			("prior",							po::value<std::string>()->default_value("prior.xml"),					"File describing the priors for the parameters/initial conditions")
			("likelihood",						po::value<std::string>()->default_value("likelihood.xml"),				"File describing the likelihood functions")
			("learning_rate,e",					po::value<Real>()->default_value(1.0),									"Learning rate between 0 and 1 - raises the likelihood function to this power.")
			("output.folder",					po::value<std::string>()->default_value("output"),						"Folder in which output files will be generated, subfolder of the model directory.")
			("predict",																									"Make a prediction for the function described in the likelihood file, given a set of output samples.")
			("predict.input",					po::value<std::string>()->default_value("output.nc"),					"This should point to the output.nc file of a previous inference.")
			("predict.output",					po::value<std::string>()->default_value("prediction.nc"),				"Prediction output will be written to this file (in the output folder).")
			("predict.skip_n",					po::value<size_t>()->default_value(0),									"Prediction output will only make output for every n+1 sample.")
			("predict.specific_temperature",	po::value<size_t>()->default_value(std::numeric_limits<size_t>::max()),	"Prediction output will only make output for every this specific temperature.")
			("progress_update_time",			po::value<Real>()->default_value(0.5),									"Update the progress counter every x seconds.")
			;
#if 0
			("predict", "Make a prediction for the function described in the likelihood file, given a set of output samples.")
			("predict.input", po::value<std::string>()->default_value("output.nc"), "This should point to the output.nc file of a previous inference.")
			("predict.output", po::value<std::string>()->default_value("prediction.nc"), "Prediction output will be written to this file (in the output folder).")
			("tempered_prediction", po::value<std::string>(), "Option for generating predictions from tempered output samples.")
			("bootstrap_celllines", po::value<unsigned long>()->default_value(0), "Do inference with a bootstrap over cell lines")
			("safe_bayes_ix", po::value<size_t>()->default_value(std::numeric_limits<size_t>::max()), "SafeBayes - include cell lines up to and including index")
			("safe_bayes_cv_ix", po::value<std::string>()->default_value(std::string()), "Cross-validated SafeBayes - include cell lines excluding index")
			;
		//bcm::SamplerFactory::AddOptionsDescription(config_options);
#endif
		bcm3::SamplerFactory::AddOptionsDescription(config_options);
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
			} else 
#endif
			if (vm.count("predict")) {
				result = predict(vm);
			} else {
				result = run(vm);
			}
		}
	} catch (po::error & e) {
		LOGERROR("Error processing program options: %s", e.what());
	}

	return result;
}
