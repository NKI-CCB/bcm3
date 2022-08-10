#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "RNG.h"
#include "UnivariateMarginal.h"

namespace bcm3 {

UnivariateMarginal::UnivariateMarginal()
	: PriorType(D_Invalid)
	, mu(0)
	, sigma(0)
	, a(0)
	, b(0)
	, lambda(0)
	, k(0)
	, theta(0)
	, scale(0)
	, lambda2(0)
	, mixing_coeff(0.5)
	, OneOverTwoSigmaSq(0)
	, LogOneOverSqrtTwoPiSigmaSq(0)
{
}

bool UnivariateMarginal::Initialize(boost::property_tree::ptree pt)
{
	try {
		std::string distribution_name = pt.get<std::string>("<xmlattr>.distribution");
		if(distribution_name == "uniform") {
			PriorType = D_Uniform;
			a = pt.get<Real>("<xmlattr>.lower");
			b = pt.get<Real>("<xmlattr>.upper");
			if (b <= a) {
				LOGERROR("Uniform distribution with upper bound less than or equal to lower bound.");
				return false;
			}
		} else if(distribution_name == "normal") {
			PriorType = D_Normal;
			mu = pt.get<Real>("<xmlattr>.mu");
			sigma = pt.get<Real>("<xmlattr>.sigma");
			OneOverTwoSigmaSq = 1.0 / (2.0 * sigma * sigma);
			LogOneOverSqrtTwoPiSigmaSq = log(rsqrt(2.0 * sigma * sigma * M_PI));
			if (sigma <= 0.0) {
				LOGERROR("Normal distribution with non-positive sigma.");
				return false;
			}
		} else if(distribution_name == "exponential") {
			PriorType = D_Exponential;
			lambda = pt.get<Real>("<xmlattr>.lambda");
			if (lambda <= 0.0) {
				LOGERROR("Exponential distribution with non-positive lambda.");
				return false;
			}
		} else if(distribution_name == "gamma") {
			PriorType = D_Gamma;
			k = pt.get<Real>("<xmlattr>.k");
			theta = pt.get<Real>("<xmlattr>.theta");
			if (k <= 0.0) {
				LOGERROR("Gamma distribution with non-positive k.");
				return false;
			}
			if (theta <= 0.0) {
				LOGERROR("Gamma distribution with non-positive theta.");
				return false;
			}
		} else if(distribution_name == "beta") {
			PriorType = D_Beta;
			a = pt.get<Real>("<xmlattr>.a");
			b = pt.get<Real>("<xmlattr>.b");
			if (a <= 0.0) {
				LOGERROR("Beta distribution with non-positive a.");
				return false;
			}
			if (b <= 0.0) {
				LOGERROR("Beta distribution with non-positive b.");
				return false;
			}
		} else if(distribution_name == "half_cauchy") {
			PriorType = D_HalfCauchy;
			scale = pt.get<Real>("<xmlattr>.scale");
			if (scale <= 0.0) {
				LOGERROR("Beta distribution with non-positive scale.");
				return false;
			}
		} else if (distribution_name == "beta_prime") {
			PriorType = D_BetaPrime;
			a = pt.get<Real>("<xmlattr>.a");
			b = pt.get<Real>("<xmlattr>.b");
			scale = pt.get<Real>("<xmlattr>.scale");
		} else if (distribution_name == "exponential_mix") {
			PriorType = D_ExponentialMix;
			lambda = pt.get<Real>("<xmlattr>.lambda");
			lambda2 = pt.get<Real>("<xmlattr>.lambda2");
			mixing_coeff = pt.get<Real>("<xmlattr>.mix");
		} else {
			PriorType = D_Invalid;
			LOGERROR("Invalid distribution type \"%s\"", distribution_name.c_str());
			return false;
		}
	} catch (boost::property_tree::ptree_error &e) {
		LOGERROR("Error parsing UnivariateMarginal: %s", e.what());
		return false;
	}

	return true;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateUniform(Real lower, Real upper)
{
	if (upper < lower) {
		LOGERROR("Incorrect uniform distribution specification for UnivariateMarginal: upper bound is smaller than the lower bound");
		return NULL;
	}

	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_Uniform;
	m->a = lower;
	m->b = upper;
	return m;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateNormal(Real mu, Real sigma)
{
	if (sigma < 0.0) {
		LOGERROR("Incorrect normal distribution specification for UnivariateMarginal: negative sigma");
		return NULL;
	}
	
	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_Normal;
	m->mu = mu;
	m->sigma = sigma;
	m->OneOverTwoSigmaSq = 1.0 / (2.0 * sigma * sigma);
	m->LogOneOverSqrtTwoPiSigmaSq = log(rsqrt(2.0 * sigma * sigma * M_PI));
	return m;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateExponential(Real lambda)
{
	if (lambda < 0.0) {
		LOGERROR("Incorrect exponential distribution specification for UnivariateMarginal: negative lambda");
		return NULL;
	}

	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_Exponential;
	m->lambda = lambda;
	return m;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateGamma(Real k, Real theta)
{
	if (k < 0.0) {
		LOGERROR("Incorrect gamma distribution specification for UnivariateMarginal: negative k");
		return NULL;
	}
	if (theta < 0.0) {
		LOGERROR("Incorrect gamma distribution specification for UnivariateMarginal: negative theta");
		return NULL;
	}

	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_Gamma;
	m->k = k;
	m->theta = theta;
	return m;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateBeta(Real a, Real b)
{
	if (a < 0.0) {
		LOGERROR("Incorrect beta distribution specification for UnivariateMarginal: negative a");
		return NULL;
	}
	if (b < 0.0) {
		LOGERROR("Incorrect beta distribution specification for UnivariateMarginal: negative b");
		return NULL;
	}
	
	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_Beta;
	m->a = a;
	m->b = b;
	return m;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateHalfCauchy(Real scale)
{
	if (scale < 0.0) {
		LOGERROR("Incorrect half-cauchy distribution specification for UnivariateMarginal: negative scale");
		return NULL;
	}
	
	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_HalfCauchy;
	m->scale = scale;
	return m;
}

std::unique_ptr<UnivariateMarginal> UnivariateMarginal::CreateBetaPrime(Real a, Real b, Real scale)
{
	if (a <= 0.0) {
		LOGERROR("Incorrect beta distribution specification for UnivariateMarginal: non-positive a");
		return NULL;
	}
	if (b <= 0.0) {
		LOGERROR("Incorrect beta distribution specification for UnivariateMarginal: non-positive b");
		return NULL;
	}
	if (scale <= 0.0) {
		LOGERROR("Incorrect half-cauchy distribution specification for UnivariateMarginal: negative scale");
		return NULL;
	}

	std::unique_ptr<UnivariateMarginal> m = std::make_unique<UnivariateMarginal>();
	m->PriorType = D_BetaPrime;
	m->a = a;
	m->b = b;
	m->scale = scale;
	return m;
}

bool UnivariateMarginal::Sample(Real& value, RNG* rng) const
{
	switch (PriorType) {
	case D_Uniform:
		value = rng->GetUniform(a, b);
		return true;

	case D_Normal:
		value = rng->GetNormal(mu, sigma);
		return true;

	case D_Exponential:
		value = rng->GetExponential(1.0 / lambda);
		return true;

	case D_Gamma:
		value = rng->GetGamma(k, theta);
		return true;

	case D_Beta:
		value = rng->GetBeta(a, b);
		return true;

	case D_HalfCauchy:
		value = rng->GetHalfCauchy(scale);
		return true;

	case D_BetaPrime:
		{
			Real x = rng->GetBeta(a, b);
			value = scale * ((x) / (1.0 - x));
		}
		return true;

	case D_ExponentialMix:
		{
			Real p = rng->GetReal();
			if (p < mixing_coeff) {
				value = rng->GetExponential(1.0 / lambda);
			} else {
				value = rng->GetExponential(1.0 / lambda2);
			}
		}
		return true;

	default:
		LOGERROR("Cannot sample from prior type %d", (int)PriorType);
		return false;
	}
}

bool UnivariateMarginal::EvaluatePDF(Real& p, const Real value) const
{
	switch (PriorType) {
	case D_Uniform:
		if (value < a || value > b) {
			p = 0.0;
		} else {
			p = 1.0 / (b - a);
		}
		return true;

	case D_Normal:
		p = PdfNormal(value, mu, sigma);
		return true;

	case D_Exponential:
		p = PdfExponential(value, lambda);
		return true;

	case D_Gamma:
		p = PdfGamma(value, k, theta);
		return true;

	case D_Beta:
		p = PdfBeta(value, a, b);
		return true;

	case D_HalfCauchy:
		if (value <= 0.0) {
			p = 0.0;
		} else {
			Real xovers = value / scale;
			p = 2.0 / (scale * M_PI * (1 + xovers * xovers));
		}
		return true;

	case D_BetaPrime:
		if (value <= 0.0) {
			p = 0.0;
		} else {
			p = PdfBetaPrime(value, a, b, scale);
		}
		return true;

	case D_ExponentialMix:
		p = mixing_coeff * PdfExponential(value, lambda) + (1.0 - mixing_coeff) * PdfExponential(value, lambda2);
		return true;

	default:
		LOGERROR("Cannot evaluate PDF for prior type %d", (int)PriorType);
		p = 0.0;
		return false;
	}
}

bool UnivariateMarginal::EvaluateLogPDF(Real& logp, const Real value) const
{
	switch (PriorType) {
	case D_Uniform:
		if (value < a || value > b) {
			logp = -std::numeric_limits<Real>::infinity();
		} else {
			logp = -log(b - a);
		}
		return true;

	case D_Normal:
		{
			const Real d = value - mu;
			logp = LogOneOverSqrtTwoPiSigmaSq - d * d * OneOverTwoSigmaSq;
		}
		return true;

	case D_Exponential:
		logp = LogPdfExponential(value, lambda);
		return true;

	case D_Gamma:
		logp = log(PdfGamma(value, k, theta));
		return true;

	case D_Beta:
		logp = log(PdfBeta(value, a, b));
		return true;

	case D_HalfCauchy:
		if (value <= 0.0) {
			logp = -std::numeric_limits<Real>::infinity();
		} else {
			logp = -0.45158270528945486472619522989488 - log(scale + value * value / scale); // -0.45158 = log(2) - log(pi)
		}
		return true;

	case D_BetaPrime:
		logp = LogPdfBetaPrime(value, a, b, scale);
		return true;

	case D_ExponentialMix:
		logp = logsum(log(mixing_coeff) + LogPdfExponential(value, lambda), log(1.0 - mixing_coeff) + LogPdfExponential(value, lambda2));
		return true;

	default:
		LOGERROR("Cannot evaluate PDF for prior type %d", (int)PriorType);
		logp = -std::numeric_limits<Real>::infinity();
		return false;
	}
}

bool UnivariateMarginal::EvaluateLogDerivative(Real& d, const Real value) const
{
	switch (PriorType) {
	case D_Uniform:
		d = 0.0;
		return true;

	case D_Normal:
		d = -2.0 * (value - mu) * OneOverTwoSigmaSq;
		return true;

	case D_Exponential:
		if (value >= 0) {
			d = -lambda;
		} else {
			d = 0;
		}
		return true;

	case D_Gamma:
	case D_Beta:
	case D_HalfCauchy:
	case D_BetaPrime:
	case D_ExponentialMix:
		LOGERROR("Derivative not implemented for prior type %d", (int)PriorType);
		d = std::numeric_limits<Real>::quiet_NaN();
		return false;

	default:
		LOGERROR("Cannot evaluate derivative for prior type %d", (int)PriorType);
		d = std::numeric_limits<Real>::quiet_NaN();
		return false;
	}
}

bool UnivariateMarginal::EvaluateLog2ndDerivative(Real& d, const Real value) const
{
	switch (PriorType) {
	case D_Uniform:
		d = 0.0;
		return true;

	case D_Normal:
		d = -2.0 * OneOverTwoSigmaSq; 
		return true;

	case D_Exponential:
		d = 0;
		return true;

	case D_Gamma:
	case D_Beta:
	case D_HalfCauchy:
	case D_BetaPrime:
	case D_ExponentialMix:
		LOGERROR("Second derivative not implemented for prior type %d", (int)PriorType);
		d = std::numeric_limits<Real>::quiet_NaN();
		return false;

	default:
		LOGERROR("Cannot evaluate gradient for prior type %d", (int)PriorType);
		d = std::numeric_limits<Real>::quiet_NaN();
		return false;
	}
}

bool UnivariateMarginal::EvaluateMean(Real& mean) const
{
	switch (PriorType) {
	case D_Uniform:
		mean = 0.5 * (b + a);
		return true;

	case D_Normal:
		mean = mu;
		return true;

	case D_Exponential:
		mean = 1.0 / lambda;
		return true;

	case D_Gamma:
		mean = k * theta;
		return true;

	case D_Beta:
		mean = a / (a + b);
		return true;

	case D_HalfCauchy:
		// Undefined
		mean = scale;
		return true;

	case D_BetaPrime:
		if (b > 1.0) {
			mean = scale * a / (b - 1);
		} else {
			// Undefined
			mean = scale;
		}
		return true;

	case D_ExponentialMix:
		mean = mixing_coeff / lambda + (1.0 - mixing_coeff) / lambda2;
		return true;

	default:
		LOGERROR("Cannot evaluate mean for prior type %d", (int)PriorType);
		return false;
	}
}

bool UnivariateMarginal::EvaluateVariance(Real& var) const
{
	switch (PriorType) {
	case D_Uniform:
		{
			Real d = b - a;
			var = (d * d) / 12.0;
		}
		return true;

	case D_Normal:
		var = sigma * sigma;
		return true;

	case D_Exponential:
		var = 1.0 / (lambda * lambda);
		return true;

	case D_Gamma:
		var = k * theta * theta;
		return true;

	case D_Beta:
		{
			Real apb = a + b;
			var = (a * b) / (apb * apb * (apb + 1));
		}
		return true;

	case D_HalfCauchy:
		// This is undefined. We'll give the squared scale as something of an approximation (this is the variance of a T-distribution with dof -> infinity).
		var = scale * scale;
		return true;

	case D_BetaPrime:
		if (b > 2.0) {
			var = scale * scale * a * (a + b - 1.0) / ((b - 2) * (b - 1) * (b - 1));
		} else {
			// Undefined
			var = scale * scale;
		}
		return true;

	case D_ExponentialMix:
		var = square(mixing_coeff) / (lambda * lambda) + square(1.0 - mixing_coeff) / (lambda2 * lambda2);
		return true;

	default:
		LOGERROR("Cannot evaluate variance for prior type %d", (int)PriorType);
		return false;
	}
}

bool UnivariateMarginal::EvaluateCDF(Real& p, const Real value) const
{
	switch (PriorType) {
	case D_Uniform:
		if (value < a || value > b) {
			p = 0.0;
		} else {
			p = (value - a) / (b - a);
		}
		return true;

	case D_Normal:
		p = CdfNormal(value, mu, sigma);
		return true;

	case D_Exponential:
		p = CdfExponential(value, lambda);
		return true;

	case D_Gamma:
		p = CdfGamma(value, k, theta);
		return true;

	case D_Beta:
		p = CdfBeta(value, a, b);
		return true;

	case D_HalfCauchy:
		p = CdfHalfCauchy(value, scale);
		return true;

	default:
		LOGERROR("Cannot evaluate cdf for prior type %d", (int)PriorType);
		p = 0.0;
		return false;
	}
}

bool UnivariateMarginal::EvaluateQuantile(Real& value, const Real p) const
{
	if (p < 0.0 || p > 1.0) {
		value = std::numeric_limits<Real>::quiet_NaN();
		return true;
	}

	switch (PriorType) {
	case D_Uniform:
		value = QuantileUniform(p, a, b);
		return true;

	case D_Normal:
		value = QuantileNormal(p, mu, sigma);
		return true;

	case D_Exponential:
		value = QuantileExponential(p, lambda);
		return true;

	case D_Gamma:
		value = QuantileGamma(p, k, theta);
		return true;

	case D_Beta:
		value = QuantileBeta(p, a, b);
		return true;

	case D_HalfCauchy:
		value = QuantileHalfCauchy(p, scale);
		return true;

	default:
		LOGERROR("Cannot evaluate quantile for prior type %d", (int)PriorType);
		value = std::numeric_limits<Real>::quiet_NaN();
		return false;
	}
}

Real UnivariateMarginal::GetLowerBound() const
{
	if (PriorType == D_Uniform) {
		return a;
	} else if (PriorType == D_Beta || PriorType == D_Exponential || PriorType == D_Gamma || PriorType == D_HalfCauchy || PriorType == D_BetaPrime) {
		return 0.0;
	} else {
		return -std::numeric_limits<Real>::infinity();
	}
}

Real UnivariateMarginal::GetUpperBound() const
{
	if (PriorType == D_Uniform) {
		return b;
	} else if (PriorType == D_Beta) {
		return 1.0;
	} else {
		return std::numeric_limits<Real>::infinity();
	}
}

}
