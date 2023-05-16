#include "Utils.h"
#include "GMM.h"
#include "NetCDFBundler.h"
#include "Prior.h"
#include "ProposalGaussianMixtureFitInR.h"
#include "SummaryStats.h"

#include <boost/process.hpp>
#include <iostream>

namespace bcm3 {

	std::mutex ProposalGaussianMixtureFitInR::netcdf_mutex;

	ProposalGaussianMixtureFitInR::ProposalGaussianMixtureFitInR()
	{
	}

	ProposalGaussianMixtureFitInR::~ProposalGaussianMixtureFitInR()
	{
	}

	bool ProposalGaussianMixtureFitInR::InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info)
	{
		if (history.rows() < 2) {
			// Start with vector sampler with a single component with diagonal covariance equal to prior variance
			VectorReal mean(num_variables);
			MatrixReal covariance(num_variables, num_variables);

			covariance.setZero();
			for (size_t i = 0; i < variable_indices.size(); i++) {
				Real prior_mean;
				if (!prior->EvaluateMarginalMean(variable_indices[i], prior_mean)) {
					// TODO - somehow come up with something reasonable?
					prior_mean = 0.0;
				}

				Real prior_var;
				if (!prior->EvaluateMarginalVariance(variable_indices[i], prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}

				mean(i) = prior_mean;
				covariance(i, i) = prior_var;
			}

			std::vector<VectorReal> means(1, mean);
			std::vector<MatrixReal> covs(1, covariance);
			VectorReal weights;
			weights.setConstant(1, 1);

			gmm = std::make_shared<GMM>();
			if (!gmm->Set(means, covs, weights)) {
				return false;
			}
		} else {
			NetCDFBundler bundler;

			LOG("Writing history to %s", tmpfilename.c_str());

			if (boost::filesystem::exists(tmpfilename)) {
				boost::filesystem::remove(tmpfilename);
			}

			{
				std::lock_guard<std::mutex> guard(netcdf_mutex);
				if (!bundler.Open(tmpfilename)) {
					LOGERROR("Unable to open file \"%s\" to store history samples temporarily", tmpfilename.c_str());
					return false;
				} else {
					std::string group = "history";
					bundler.AddGroup(group);
					bundler.AddMatrix(group, "samples", history);
					bundler.Close();
				}
			}

			std::string bcm_path;
			char* bcm_root_env = getenv("BCM3_ROOT");
			if (!bcm_root_env) {
				LOGERROR("BCM3_ROOT environment variable has not been specified!");
				return false;
			} else {
				bcm_path = bcm_root_env;
				bcm3::fix_path(bcm_path);
			}

			std::string cmd = std::string("Rscript ") + bcm_path + std::string("R/fit_proposal.r ") + tmpfilename.c_str();
			if (log_info) {
				cmd += std::string(" log_info");
			}
			LOG("Calling R fitting script using: \"%s\"", cmd.c_str());
			int result = boost::process::system(cmd);

			if (result != 0) {
				LOGERROR("R script for fitting proposal distribution failed");
				return false;
			}

			{
				std::lock_guard<std::mutex> guard(netcdf_mutex);

				NetCDFDataFile fit_output;
				if (!fit_output.Open(tmpfilename, false)) {
					LOGERROR("Unable to open file \"%s\" to read fitted proposal", tmpfilename.c_str());
					return false;
				} else {
					VectorReal weights;
					fit_output.GetVector("fitted_proposal", "weights", weights);

					std::vector<VectorReal> means(weights.size());
					std::vector<MatrixReal> covariances(weights.size());
					for (ptrdiff_t i = 0; i < weights.size(); i++) {
						fit_output.GetVector("fitted_proposal", std::string("mean") + std::to_string(i + 1), means[i]);
						fit_output.GetMatrix("fitted_proposal", std::string("covariance") + std::to_string(i + 1), covariances[i]);
					}

					gmm = std::make_shared<GMM>();
					gmm->Set(means, covariances, weights);

					fit_output.Close();
					boost::filesystem::remove(tmpfilename);
				}
			}
		}

		scales.setConstant(gmm->GetNumComponents(), 2.38 / sqrt(num_variables));
		acceptance_rate_emas.setConstant(gmm->GetNumComponents(), target_acceptance_rate);

		return true;
	}

}
