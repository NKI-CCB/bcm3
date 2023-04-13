#include "Utils.h"
#include "SamplerPT.h"
#include "VariableSet.h"

namespace bcm3 {

SamplerPT::SamplerPT(size_t threads, size_t max_memory_use)
	: numthreads(threads)
	, num_samples(25000)
	, use_every_nth(1)
	, num_variables(0)
	, blocking_strategy("one_block")
	, proposal_type("global_covariance")
	, rngseed(0)
	, output_proposal_adaptation(false)
	, swapping_scheme(ESwappingScheme::StochasticEvenOdd)
	, exchange_probability(0.5)
	, num_exploration_steps(1)
	, adapt_proposal_samples(5000)
	, adapt_proposal_times(2)
	, history_max_samples(20000)
	, adapt_proposal_max_samples(2000)
	, proposal_scaling_learning_rate(0.05)
	, proposal_scaling_ema_period(5000)
	, stop_proposal_scaling(15000)
	, proposal_adaptations_done(0)
	, proposal_scaling_adaptations_done(false)
	, previous_swap_even(false)
	, target_acceptance_rate(0.23)
	, proposal_t_df(0.0)
	, num_likelihood_evaluations(0)
{
	parallel_likelihoods.resize(threads);
	size_t default_num_chains = 6;
	fixed_temperatures.setZero(default_num_chains);
	for (size_t i = 1; i < default_num_chains; i++) {
		Real alpha = i / (Real)(default_num_chains - 1);
		fixed_temperatures(i) = pow(alpha, 3.0);
	}
	fixed_temperatures(default_num_chains-1) = 1.0;
}

SamplerPT::~SamplerPT()
{
}

void SamplerPT::SetVariableSet(std::shared_ptr<VariableSet> varset)
{
	this->varset = varset;
	num_variables = varset->GetNumVariables();
}

void SamplerPT::SetPrior(std::shared_ptr<Prior> prior)
{
	this->prior = prior;
}

void SamplerPT::SetLikelihood(std::shared_ptr<Likelihood> likelihood, size_t threadix)
{
	parallel_likelihoods[threadix] = likelihood;
}

void SamplerPT::SetOutputPath(const std::string& path)
{
	output_path = path;
}

void SamplerPT::SetProgressIndicator(std::shared_ptr<ProgressIndicator> pi)
{
	this->progress = pi;
}

bool SamplerPT::LoadSettings(const boost::program_options::variables_map& vm)
{
	try {
		num_samples = vm["sampler.num_samples"].as<size_t>();
		use_every_nth = vm["sampler.use_every_nth"].as<size_t>();
		blocking_strategy = vm["sampler.blocking_strategy"].as<std::string>();
		proposal_type = vm["sampler.proposal_type"].as<std::string>();
		history_max_samples = vm["sampler.history_max_samples"].as<size_t>();
		adapt_proposal_samples = vm["sampler.adapt_proposal_samples"].as<size_t>();
		adapt_proposal_times = vm["sampler.adapt_proposal_times"].as<size_t>();
		adapt_proposal_max_samples = vm["sampler.adapt_proposal_max_samples"].as<size_t>();
		exchange_probability = vm["sampler.exchange_probability"].as<Real>();
		num_exploration_steps = vm["sampler.num_exploration_steps"].as<size_t>();
		proposal_t_df = vm["sampler.proposal_t_df"].as<Real>();
		clustered_blocking_nn = vm["sampler.clustered_blocking_nn"].as<int>();
		clustered_blocking_nn2 = vm["sampler.clustered_blocking_nn2"].as<int>();
		clustered_blocking_n_clusters = vm["sampler.clustered_blocking_n_clusters"].as<int>();
		rngseed = vm["sampler.rngseed"].as<unsigned long long>();
		output_proposal_adaptation = vm["sampler.output_proposal_adaptation"].as<bool>();

		std::string swapping_scheme_str = vm["sampler.swapping_scheme"].as<std::string>();
		if (swapping_scheme_str == "stochastic_random") {
			swapping_scheme = ESwappingScheme::StochasticRandom;
		} else if (swapping_scheme_str == "stochastic_even_odd") {
			swapping_scheme = ESwappingScheme::StochasticEvenOdd;
		} else if (swapping_scheme_str == "deterministic_even_odd") {
			swapping_scheme = ESwappingScheme::DeterministicEvenOdd;
		} else {
			LOGERROR("Unknown swapping scheme \"%s\"", swapping_scheme_str.c_str());
			return false;
		}

		if (adapt_proposal_samples > 0) {
			proposal_scaling_ema_period = (size_t)ceil(adapt_proposal_samples * use_every_nth * (1.0 - exchange_probability) / 10);
			proposal_scaling_learning_rate = pow(100.0, 1.0 / proposal_scaling_ema_period) - 1.0;
		}

		size_t num_chains = vm["sampler.num_chains"].as<size_t>();
		Real temperature_power = vm["sampler.fixed_temperature_schedule_power"].as<Real>();
		Real temperature_max = vm["sampler.fixed_temperature_schedule_max"].as<Real>();

		// Power-law scaling from 1.0 downwards and the lowest chain at temperature 0
		fixed_temperatures.setZero(num_chains);
		for (size_t i = 1; i < num_chains - 1; i++) {
			Real alpha = i / (Real)(num_chains - 1);
			fixed_temperatures(i) = temperature_max * pow(alpha, temperature_power);
		}
		fixed_temperatures(num_chains - 1) = temperature_max;
	} catch (boost::program_options::error& e) {
		LOGERROR("Error parsing sampler: %s", e.what());
		return false;
	}

	return true;
}

bool SamplerPT::Initialize()
{
	if (fixed_temperatures.size() == 0) {
		LOGERROR("Temperatures have not been specified.");
		return false;
	}
	if (!varset) {
		LOGERROR("No variable set specified");
		return false;
	}
	if (num_variables == 0) {
		LOGERROR("No variables to sample");
		return false;
	}

	// Start sampling with chains at the fixed temperatures
	current_temperatures = fixed_temperatures;

	// Calculate the number of samples we can store in history
	size_t expected_sample_history = adapt_proposal_samples * use_every_nth;
	size_t history_subsampling = 1;
	size_t sample_history = expected_sample_history;
	if (sample_history > history_max_samples) {
		history_subsampling = (expected_sample_history + history_max_samples - 1) / history_max_samples;
		sample_history = adapt_proposal_samples * use_every_nth / history_subsampling;
	}
	if (adapt_proposal_times > 0) {
		LOG("Using history size %zd with 1 in %zd subsampling for proposal adaptation", sample_history, history_subsampling);
	}

	// Initialize chains
	chains.resize(current_temperatures.size());
	for (size_t i = 0; i < chains.size(); i++) {
		chains[i] = std::make_unique<SamplerPTChain>(this);
		chains[i]->SetTemperature(current_temperatures(i), true);
		bool result = chains[i]->Initialize(sample_history, history_subsampling);
		if (!result) {
			return false;
		}
	}
	previous_swap_even = false;

	// Initialize proposal adaptations
	proposal_adaptations_done = 0;
	proposal_scaling_adaptations_done = false;

	// Initialize multi-threading structure
	task_manager = std::make_unique<TaskManager>(numthreads);
	async.resize(numthreads);

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

void SamplerPT::AddSampleHandler(std::shared_ptr<SampleHandler> handler)
{
	sample_handlers.push_back(handler);
}

bool SamplerPT::Run()
{
	bool result = true;
	Real max_lposterior = -std::numeric_limits<Real>::infinity();
	num_likelihood_evaluations = 0;

	bcm3::Timer timer;
	timer.Start();
	if (auto pr = progress.lock()) {
		pr->NotifyStart();
	}

	LOG("Finding starting point...");
	for (auto& chain : chains) {
		result = chain->FindStartingPosition();
		if (!result) {
			return false;
		}
	}
	
	LOG("Starting sampling loop...");
	UpdateProgress(0.0, true);
	size_t total_samples = num_samples * use_every_nth;
	for (size_t si = 0; si < total_samples; si++) {
		UpdateProgress(si / (Real)total_samples, false);
		if (total_samples > 100) {
			if (si % (total_samples / 100) == 0) {
				double progress = 100.0 * si / (Real)total_samples;
				LOG("Progress %.0f%%", progress);
			}
		} else {
			double progress = 100.0 * si / (Real)total_samples;
			LOG("Progress %.0f%%", progress);
		}

		size_t sample_ix = si / use_every_nth;

		if (chains.size() > 1) {
			if (swapping_scheme == ESwappingScheme::StochasticRandom || swapping_scheme == ESwappingScheme::StochasticEvenOdd) {
				// Randomly choose between mutate and exchange moves
				Real exchange = rng.GetReal();
				if (exchange < exchange_probability) {
					DoExchangeMove(si);
				} else {
					result = DoMutateMove();
				}
			} else if (swapping_scheme == ESwappingScheme::DeterministicEvenOdd) {
				for (size_t ei = 0; ei < num_exploration_steps; ei++) {
					result = DoMutateMove();
					if (!result) {
						break;
					}
				}
				DoExchangeMove(si);
			}
		} else {
			result = DoMutateMove();
		}

		if (!result) {
			LOGERROR("Sample step failed");
			LogStatistics();
			return false;
		}

		if ((*chains.rbegin())->GetLogPowerPosterior() > max_lposterior) {
			max_lposterior = (*chains.rbegin())->GetLogPowerPosterior();
			LOG("New max log posterior: %g", max_lposterior);
		}

		if ((si + 1) % use_every_nth == 0) {
			EmitSample(sample_ix);

			if (adapt_proposal_samples > 0 && ((sample_ix + 1) % adapt_proposal_samples == 0) && si + 1 != total_samples && proposal_adaptations_done < adapt_proposal_times) {
				// Force the progress update
				UpdateProgress(si / (Real)total_samples, true);

				LogStatistics();
				LOG("Updating proposal to empirical variance");
				bool result = true;
				for (auto& chain : chains) {
					chain->AdaptProposalAsync();
				}
				for (auto& chain : chains) {
					result &= chain->AdaptProposalWait();
				}
				proposal_adaptations_done++;
			}
			if (stop_proposal_scaling > 0 && sample_ix > stop_proposal_scaling) {
				if (!proposal_scaling_adaptations_done) {
					LOG("Stopping adaptation of proposal scale");
					proposal_scaling_adaptations_done = true;
				}
			}
		}
	}

	UpdateProgress(1.0, true);
	if (auto pr = progress.lock()) {
		pr->NotifyStop();
	}

	LOG("Finished!");
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

void SamplerPT::AddOptionsDescription(boost::program_options::options_description& pod)
{
	pod.add_options()
		("sampler.num_samples",							boost::program_options::value<size_t>()->default_value(2500),							"Number of samples to produce, after the burnin period.")
		("sampler.use_every_nth",						boost::program_options::value<size_t>()->default_value(1),								"Subsampling, a value of 10 means that every 10th sample is used and the 9 samples in between are discarded.")
		("sampler.num_chains",							boost::program_options::value<size_t>()->default_value(6),								"Number of fixed temperature chains")
		("sampler.blocking_strategy",					boost::program_options::value<std::string>()->default_value("autoblock"),				"Blocking strateg: one_block, no_blocking, Turek")
		("sampler.proposal_type",						boost::program_options::value<std::string>()->default_value("global_covariance"),		"Type of proposal: global_covariance")
		("sampler.adapt_proposal_samples",				boost::program_options::value<size_t>()->default_value(500),							"Number of samples after which the proposal distribution should be adapted, 0 for no adaptation.")
		("sampler.adapt_proposal_times",				boost::program_options::value<size_t>()->default_value(2),								"Number of times the proposal variance should be adapted.")
		("sampler.history_max_samples",					boost::program_options::value<size_t>()->default_value(20000),							"Maximum number of samples to store in the sample history.")
		("sampler.adapt_proposal_max_samples",			boost::program_options::value<size_t>()->default_value(2000),							"Maximum number of samples to use in the GMM fitting/spectral clustering.")
		("sampler.stop_proposal_scaling",				boost::program_options::value<size_t>()->default_value(1500),							"Stop adaptive scaling of the proposal distribution after this many samples.")
		("sampler.swapping_scheme",						boost::program_options::value<std::string>()->default_value("deterministic_even_odd"),	"Swapping scheme, can be either stochastic_random, stochastic_even_odd or deterministic_even_odd.")
		("sampler.exchange_probability",				boost::program_options::value<Real>()->default_value(0.5),								"Probability of choosing an exchange move instead of a mutate move (only used if swapping_scheme is stochastic).")
		("sampler.num_exploration_steps",				boost::program_options::value<size_t>()->default_value(1),								"Number of exploration steps between swaps (only used if swapping_scheme is deterministic).")
		("sampler.proposal_t_df",						boost::program_options::value<Real>()->default_value(0.0),								"Use a t-distribution with the specified degrees of freedom as proposal (0 for a normal distribution).")
		("sampler.clustered_blocking_nn",				boost::program_options::value<int>()->default_value(15),								"")
		("sampler.clustered_blocking_nn2",				boost::program_options::value<int>()->default_value(30),								"")
		("sampler.clustered_blocking_n_clusters",		boost::program_options::value<int>()->default_value(4),									"")
		("sampler.fixed_temperature_schedule_power",	boost::program_options::value<Real>()->default_value(3.0),								"Specifies the rate of the power-law used for the initial temperature schedule.")
		("sampler.fixed_temperature_schedule_max",		boost::program_options::value<Real>()->default_value(1.0),								"Specifies the maximum temperature used for the initial temperature schedule.")
		("sampler.rngseed",								boost::program_options::value<unsigned long long>()->default_value(0),					"Random number generator seed, set to 0 to use an arbitrary seed based on the computer time.")
		("sampler.output_proposal_adaptation",			boost::program_options::value<bool>()->default_value(false),							"Whether to output information describing the proposal adaptation.")
		;
}

bool SamplerPT::EvaluatePrior(size_t threadix, const VectorReal& values, Real& lprior)
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

bool SamplerPT::EvaluateLikelihood(size_t threadix, const VectorReal& values, Real& llh)
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

bool SamplerPT::EvaluatePriorLikelihood(size_t threadix, const VectorReal& values, Real& lprior, Real& llh)
{
	if (!EvaluatePrior(threadix, values, lprior)) {
		return false;
	}
	return EvaluateLikelihood(threadix, values, llh);
}

void SamplerPT::DoExchangeMove(size_t sample_ix)
{
	if (swapping_scheme == ESwappingScheme::StochasticEvenOdd || swapping_scheme == ESwappingScheme::DeterministicEvenOdd) {
		// Alternatingly swap the odd and even chains with the subsequent chains
		size_t start_ix;
		if (previous_swap_even) {
			start_ix = 1;
			previous_swap_even = false;
		} else {
			start_ix = 0;
			previous_swap_even = true;
		}

		size_t ci = start_ix;
		for (; ci < chains.size(); ci += 2) {
			SamplerPTChain& chain1 = *chains[ci];
			size_t ix2 = ci + 1;
			if (ix2 == chains.size()) {
				ix2 = 0;
			}
			SamplerPTChain& chain2 = *chains[ix2];
			bool exchanged = chain1.ExchangeMove(chain2);
		}
	} else if (swapping_scheme == ESwappingScheme::StochasticRandom) {
		size_t ci = rng.GetUnsignedInt(chains.size()-2);
		SamplerPTChain& chain1 = *chains[ci];
		SamplerPTChain& chain2 = *chains[ci+1];
		chain1.ExchangeMove(chain2);
	}
}

bool SamplerPT::DoMutateMove()
{
	bool result = true;
	for (auto& chain : chains) {
		chain->MutateMoveAsync();
	}
	for (auto& chain : chains) {
		result &= chain->MutateMoveWait();
	}

	return result;
}

void SamplerPT::EmitSample(size_t sample_ix)
{
	for (auto handler : sample_handlers) {
		for (auto& chain : chains) {
			if (chain->GetIsFixedTemperature()) {
				handler->ReceiveSample(chain->GetVarValues(), chain->GetLogPrior(), chain->GetLogLikelihood(), chain->GetTemperature());
			}
		}
	}
}

void SamplerPT::LogStatistics()
{
	// Statistics
	LOG("");
	LOG("Acceptance statistics:");
	LOG("Temperature | Mutate (all) | Exchange (all)");
	for (const auto& chain : chains) {
		chain->LogStatistics();
	}
	
	LOG("");
	LOG("Proposal info:");
	(*chains.rbegin())->LogProposalInfo();
}

void SamplerPT::UpdateProgress(Real p, bool force)
{
	if (auto pr = progress.lock()) {
		pr->ReportProgress(p, force);
	}
}

}
