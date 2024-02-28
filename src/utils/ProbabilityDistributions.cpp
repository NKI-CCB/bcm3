#include "Utils.h"
#include "ProbabilityDistributions.h"
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

namespace bcm3 {

Real PdfBeta(Real x, Real a, Real b)
{
	if (x < 0.0 || x > 1.0 ) {
		return 0.0;
	}
	boost::math::beta_distribution<Real> dist(a, b);
	return boost::math::pdf(dist, x);
}

Real PdfBetaPrime(Real x, Real a, Real b, Real scale)
{
	if (x < 0.0) {
		return 0.0;
	}
	Real nc = 1.0 / (boost::math::beta(a, b) * scale);
	Real sx = x / scale;
	Real rv = pow(sx, a - 1.0) * pow(sx + 1.0, -(a+b));
	return nc * rv;
}

Real PdfExponential(Real x, Real lambda)
{
	if (x < 0.0) {
		return 0.0;
	}
	return lambda * exp(-lambda * x);
}

Real PdfGamma(Real x, Real shape, Real scale)
{
	if (x < 0.0 || x == std::numeric_limits<Real>::infinity()) {
		return 0.0;
	}
	if (shape == 1.0) {
		return PdfExponential(x, 1.0 / scale);
	}
	boost::math::gamma_distribution<Real> dist(shape, scale);
	return boost::math::pdf(dist, x);
}

Real PdfNormal(Real x, Real mu, Real sigma)
{
	Real twoSigmaSq = 2.0 * sigma * sigma;
	Real d = x - mu;
	return bcm3::rsqrt(twoSigmaSq * M_PI) * exp(-(d * d) / twoSigmaSq);
}

Real PdfT(Real x, Real mu, Real sigma, Real nu)
{
	Real invsigma = 1.0 / sigma;
	Real xn = (x - mu) * invsigma;
	boost::math::students_t_distribution<Real> dist(nu);
	return boost::math::pdf(dist, xn) * invsigma;
}

Real PdfT(Real x, Real mu, Real sigma, Real nu, Real C)
{
	if (nu > 1e10) {
		return PdfNormal(x, mu, sigma);
	}
	
	Real invsigma = 1.0 / sigma;
	Real xn = (x - mu) * invsigma;
	Real basem1 = xn * xn / nu;
	Real result;
	if (basem1 < 0.125) {
		result = exp(-boost::math::log1p(basem1) * (1.0 + nu) / 2.0);
	} else {
		result = pow(1.0 / (1.0 + basem1), (nu + 1.0) / 2);
	}
	// C = 1.0 / (sqrt(nu) * boost::math::beta(nu / 2.0, 0.5))
	return invsigma * C * result;
}

Real PdfGPD(Real x, Real u, Real xi, Real beta)
{
	Real z = (x - u) / beta;
	Real p = pow(1.0 + xi * z, -1.0 / xi - 1.0) / beta;
	return p;
}

Real LogPdfBetaPrecision(Real x, Real mu, Real phi)
{
	if (mu <= 0.0 || mu >= 1.0) {
		return -std::numeric_limits<Real>::infinity();
	}
	Real a = mu * phi;
	Real b = phi * (1.0 - mu);
	if (x < 0.0 || x > 1.0) {
		return -std::numeric_limits<Real>::infinity();
	} else if (b < 1.0 && x == 1.0) {
		return -std::numeric_limits<Real>::infinity();
	} else if (a < 1.0 && x == 0.0) {
		return -std::numeric_limits<Real>::infinity();
	}
	boost::math::beta_distribution<Real> dist(a, b);
	return log(boost::math::pdf(dist, x));
}

Real LogPdfBetaPrime(Real x, Real a, Real b, Real scale)
{
	if (x < 0.0) {
		return -std::numeric_limits<Real>::infinity();
	}
	Real lnc = -log(boost::math::beta(a, b) * scale);
	Real sx = x / scale;
	Real lrv = (a - 1.0) * log(sx) - (a + b) * log(sx + 1.0);
	return lnc + lrv;
}

Real LogPdfExponential(Real x, Real lambda)
{
	if (x < 0.0) {
		return -std::numeric_limits<Real>::infinity();
	}
	return log(lambda) - lambda * x;
}

Real LogPdfNormal(Real x, Real mu, Real sigma, bool skip_na)
{
	if (skip_na && x != x) {
		return 0;
	}
	Real twoSigmaSq = 2.0 * sigma * sigma;
	Real d = x - mu;
	//return log(bcm3::rsqrt(twoSigmaSq * M_PI)) - d * d / twoSigmaSq;
	return -log(sigma) - 0.91893853320467274178032973640562 - d * d / twoSigmaSq; // 0.9189.. = log(1.0 / sqrt(2 * pi))
}

Real LogPdfTruncatedNormal(Real x, Real mu, Real sigma, Real a, Real b, bool skip_na)
{
	if (skip_na && x != x) {
		return 0.0;
	}
	if (x < a || x > b) {
		return -std::numeric_limits<Real>::infinity();
	}
	Real scale = CdfNormal(b, mu, sigma) - CdfNormal(a, mu, sigma);
	if (scale / sigma <= std::numeric_limits<Real>::epsilon()) {
		return -std::numeric_limits<Real>::infinity();
	}

	Real twoSigmaSq = 2.0 * sigma * sigma;
	Real d = x - mu;
	//return log(bcm3::rsqrt(twoSigmaSq * M_PI)) - d * d / twoSigmaSq;
	return -log(sigma * scale) - 0.91893853320467274178032973640562 - d * d / twoSigmaSq; // 0.9189.. = log(1.0 / sqrt(2 * pi))
}

Real LogPdfT(Real x, Real mu, Real sigma, Real nu, bool skip_na)
{
#if 0
	return log(PdfT(x, mu, sigma, nu));
#else
	if (nu > 1e10) {
		return LogPdfNormal(x, mu, sigma);
	}
	if (skip_na && x != x) {
		return 0;
	}
	
	Real xn = (x - mu) * sigma;
	Real basem1 = xn * xn / nu;
	if (basem1 == std::numeric_limits<Real>::infinity()) {
		return -std::numeric_limits<Real>::infinity();
	}
	Real result = -0.5 * (nu + 1.0) * boost::math::log1p(basem1);
	Real logC = -log(sigma * sqrt(nu) * boost::math::beta(nu / 2.0, 0.5));
	return logC + result;
#endif
}

Real LogPdfT(Real x, Real mu, Real sigma, Real nu, Real logC, bool skip_na)
{
	if (nu > 1e10) {
		return LogPdfNormal(x, mu, sigma);
	}
	if (skip_na && x != x) {
		return 0;
	}
	
	Real xn = (x - mu) / sigma;
	Real basem1 = xn * xn / nu;
	if (basem1 == std::numeric_limits<Real>::infinity()) {
		return -std::numeric_limits<Real>::infinity();
	}
	Real result = -0.5 * (nu + 1.0) * boost::math::log1p(basem1);
	return logC + result - log(sigma);
}

Real LogPdfT_CalcC(Real nu)
{
	Real logC = -log(sqrt(nu) * boost::math::beta(nu / 2.0, 0.5));
	return logC;
}

Real LogPdfTnu3(Real x, Real mu, Real sigma, bool skip_na)
{
	if (skip_na && x != x) {
		return 0.0;
	}

	Real xn = (x - mu) / sigma;
	return -1.0008888496235098 - 2.0 * boost::math::log1p(0.333333333333333333 * xn * xn) - log(sigma);
}

inline Real CdfTnu3(Real t)
{
	return 0.091888149236965 * ((6.0 * t) / (t * t + 3.0) + 3.464101615137754 * atan(0.577350269189626 * t) + 5.441398092702653);
}

Real LogTruncatedPdfTnu3(Real x, Real mu, Real sigma, Real a, Real b, bool skip_na)
{
	if (skip_na && x != x) {
		return 0.0;
	}
	if (x < a || x > b) {
		return -std::numeric_limits<Real>::infinity();
	}
	Real invsigma = 1 / sigma;
	Real scale = CdfTnu3((b - mu) * invsigma) - CdfTnu3((a - mu) * invsigma);
	//return LogPdfTnu3(x, mu, sigma) - log(scale);
	Real xn = (x - mu) * invsigma;
	return -1.0008888496235098 - 2.0 * boost::math::log1p(0.333333333333333333 * xn * xn) - log(sigma * scale);
}

Real LogPdfGPD(Real x, Real u, Real xi, Real beta)
{
	Real z = (x - u) / beta;
	Real logp = log(1.0 / beta) - (1.0 / xi + 1.0) * log(1.0 + xi * z);
	return logp;
}

Real CdfBeta(Real x, Real a, Real b)
{
	boost::math::beta_distribution<Real> dist(a, b);
	return boost::math::cdf(dist, x);
}

Real CdfCauchy(Real x, Real scale)
{
	return 0.31830988618379067153776752674503 * atan(x / scale) + 0.5;
}

Real CdfExponential(Real x, Real lambda)
{
	if (x < 0.0) {
		return 0.0;
	}
	return 1.0 - exp(-lambda * x);
}

Real CdfGamma(Real x, Real shape, Real scale)
{
	if (x < 0.0) {
		return 0.0;
	}
	boost::math::gamma_distribution<Real> dist(shape, scale);
	return boost::math::cdf(dist, x);
}

Real CdfHalfCauchy(Real x, Real scale)
{
	if (x < 0.0) {
		return 0.0;
	}
	return 0.63661977236758134307553505349006 * atan(x / scale);
}

Real CdfNormal(Real x, Real mu, Real sigma)
{
	boost::math::normal_distribution<Real> dist(mu, sigma);
	return boost::math::cdf(dist, x);
}

Real CdfT(Real x, Real mu, Real sigma, Real nu)
{
	boost::math::students_t_distribution<Real> dist(nu);
	return boost::math::cdf(dist, (x - mu) / sigma);
}

Real CdfGPD(Real x, Real u, Real xi, Real beta)
{
	if (x < u) {
		return 0.0;
	}
	if (xi < 0 && x > u - beta / xi) {
		return 1.0;
	}

	Real z = (x - u) / beta;
	return (1.0 - pow(1.0 + xi * z, -1.0 / xi));
}

Real CdfTruncatedNormal(Real x, Real mu, Real sigma, Real a, Real b)
{
	if (x < a || x > b) {
		return 0.0;
	}
	boost::math::normal_distribution<Real> dist(mu, sigma);
	Real cdfa = boost::math::cdf(dist, a);
	Real scale = boost::math::cdf(dist, b) - cdfa;
	return (boost::math::cdf(dist, x) - cdfa) / scale;
}

Real QuantileBeta(Real p, Real a, Real b)
{
	boost::math::beta_distribution<Real> dist(a, b);
	return boost::math::quantile(dist, p);
}

Real QuantileHalfCauchy(Real p, Real scale)
{
	if (p == 0.0) {
		return 0.0;
	} else if (p == 1.0) {
		return std::numeric_limits<Real>::infinity();
	} else {
		return scale * tan(1.5707963267948966192313216916398 * p);
	}
}

Real QuantileExponential(Real p, Real lambda)
{
	if (p == 0.0) {
		return 0.0;
	} else if (p == 1.0) {
		return std::numeric_limits<Real>::infinity();
	} else {
		return -boost::math::log1p(-p) / lambda;
	}
}

Real QuantileGamma(Real p, Real shape, Real scale)
{
	boost::math::gamma_distribution<Real> dist(shape, scale);
	return boost::math::quantile(dist, p);
}

Real QuantileNormal(Real p, Real mu, Real sigma)
{
	boost::math::normal_distribution<Real> dist(mu, sigma);
	return boost::math::quantile(dist, p);
}

Real QuantileT(Real p, Real mu, Real sigma, Real nu)
{
	boost::math::students_t_distribution<Real> dist(nu);
	return boost::math::quantile(dist, p) * sigma + mu;
}

Real QuantileUniform(Real p, Real a, Real b)
{
	if (p == 0.0) {
		return a;
	} else if (p == 1.0) {
		return b;
	} else {
		return a + p * (b - a);
	}
}

Real QuantileGPD(Real x, Real u, Real xi, Real beta)
{
	return u + beta * (pow(1.0 - x, -xi) - 1.0) / xi;
}

}
