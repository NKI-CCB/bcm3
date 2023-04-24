#include "Utils.h"
#include "Sampler.h"
#include "VariableSet.h"

namespace bcm3 {

	Sampler::Sampler(size_t threads, size_t max_memory_use)
		: num_threads(threads)
		, num_samples(25000)
		, use_every_nth(1)
		, num_variables(0)
		, logged_progress(0)
		, rngseed(0)
		, num_likelihood_evaluations(0)
	{
		parallel_likelihoods.resize(threads);
	}

	Sampler::~Sampler()
	{
	}

	void Sampler::SetVariableSet(std::shared_ptr<VariableSet> varset)
	{
		this->varset = varset;
		num_variables = varset->GetNumVariables();
	}

	void Sampler::SetPrior(std::shared_ptr<Prior> prior)
	{
		this->prior = prior;
	}

	void Sampler::SetLikelihood(std::shared_ptr<Likelihood> likelihood, size_t threadix)
	{
		parallel_likelihoods[threadix] = likelihood;
	}

	void Sampler::SetOutputPath(const std::string& path)
	{
		output_path = path;
	}

	void Sampler::SetProgressIndicator(std::shared_ptr<ProgressIndicator> pi)
	{
		progress = pi;
	}

	void Sampler::AddSampleHandler(std::shared_ptr<SampleHandler> handler)
	{
		sample_handlers.push_back(handler);
	}

	bool Sampler::LoadSettings(const boost::program_options::variables_map& vm)
	{
		try {
			num_samples = vm["sampler.num_samples"].as<size_t>();
			use_every_nth = vm["sampler.use_every_nth"].as<size_t>();
			rngseed = vm["sampler.rngseed"].as<unsigned long long>();
		} catch (boost::program_options::error& e) {
			LOGERROR("Error parsing sampler: %s", e.what());
			return false;
		}

		return true;
	}

	bool Sampler::Initialize()
	{
		if (!varset) {
			LOGERROR("No variable set specified");
			return false;
		}
		if (num_variables == 0) {
			LOGERROR("No variables to sample");
			return false;
		}

		// Initialize multi-threading structure
		task_manager = std::make_unique<TaskManager>(num_threads);
		async.resize(num_threads);

		// Initialize random number generators
		if (rngseed == 0) {
			rngseed = rng.GetTimeBasedSeed();
		}
		rng.Seed(rngseed);
		for (auto ai : async) {
			unsigned long long newseed = rng.GetUnsignedInt();
			ai.rng.Seed(newseed);
		}

		return true;
	}

	bool Sampler::Run()
	{
		bool result = true;
		num_likelihood_evaluations = 0;

		bcm3::Timer timer;
		timer.Start();
		if (auto pr = progress.lock()) {
			pr->NotifyStart();
		}
		logged_progress = 0;
		UpdateProgress(0.0, true);

		if (!RunImpl()) {
			LOGERROR("Sampling failed");
		}

		UpdateProgress(1.0, true);
		if (auto pr = progress.lock()) {
			pr->NotifyStop();
		}

		LOG("Sampling finished.");

		LogStatistics();

		double integration_time = timer.GetElapsedSeconds();
		LOG("");
		LOG("Timing information:");
		LOG("  Integration time: %.3g seconds", integration_time);
		LOG("  Number of likelihood evaluations: %u", (size_t)num_likelihood_evaluations);
		LOG("  Evaluations per second: %.0f", num_likelihood_evaluations / integration_time);
		LOG("");
		LOG("CPU information:");
		bcm3::logger->LogCPUInfo();

		return true;
	}

	void Sampler::AddOptionsDescription(boost::program_options::options_description& pod)
	{
		pod.add_options()
			("sampler.num_samples",		boost::program_options::value<size_t>()->default_value(2500),			"Number of samples to produce.")
			("sampler.use_every_nth",	boost::program_options::value<size_t>()->default_value(1),				"Subsampling, a value of 10 means that every 10th sample is used and the 9 samples in between are discarded.")
			("sampler.rngseed",			boost::program_options::value<unsigned long long>()->default_value(0),	"Random number generator seed, set to 0 to use an arbitrary seed based on the computer time. A particular seed only returns the same values when the sampling is also done with the same number of threads.")
			;
	}

	bool Sampler::EvaluatePrior(size_t threadix, const VectorReal& values, Real& lprior)
	{
		lprior = -std::numeric_limits<Real>::infinity();
		bool result = prior->EvaluateLogPDF(threadix, values, lprior);
		if (!result) {
			LOGERROR("Prior evaluation failed");
		} else if (lprior != lprior) {
			LOGERROR("NAN in prior calculation");
			result = false;
		}
		return result;
	}

	bool Sampler::EvaluateLikelihood(size_t threadix, const VectorReal& values, Real& llh)
	{
		llh = -std::numeric_limits<Real>::infinity();
		bool result = parallel_likelihoods[threadix]->EvaluateLogProbability(threadix, values, llh);
		llh *= parallel_likelihoods[threadix]->GetLearningRate();
		num_likelihood_evaluations++;
		if (!result) {
			LOGERROR("Likelihood evaluation failed");
		} else if (llh != llh) {
			LOGERROR("NAN in likelihood calculation with values:");
			for (size_t i = 0; i < num_variables; i++) {
				LOGERROR("  %s: %g", varset->GetVariableName(i).c_str(), values(i));
			}
			result = false;
		}
		return result;
	}

	bool Sampler::EvaluatePriorLikelihood(size_t threadix, const VectorReal& values, Real& lprior, Real& llh)
	{
		if (!EvaluatePrior(threadix, values, lprior)) {
			return false;
		}
		return EvaluateLikelihood(threadix, values, llh);
	}

	void Sampler::UpdateProgress(Real p, bool force)
	{
		if (auto pr = progress.lock()) {
			pr->ReportProgress(p, force);
		}

		int integer_progress = (int)(p * 100);
		if (integer_progress > logged_progress) {
			LOG("Progress %d%%", integer_progress);
			logged_progress = integer_progress;
		}
	}
}
