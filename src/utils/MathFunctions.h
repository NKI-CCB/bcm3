#pragma once

#include <boost/math/special_functions/log1p.hpp>
#include <xmmintrin.h>

namespace bcm3 {

inline Real square(Real x) { return x*x; }
inline float log2(float x) { return 1.4426950408889634073599246810019f * logf(x); }
inline double log2(double x) { return 1.4426950408889634073599246810019 * log(x); }
inline float log10(float x) { return 0.4342944819032518276511289189166f * logf(x); }
inline double log10(double x) { return 0.4342944819032518276511289189166 * log(x); }
inline float fastpow10(float x) { return exp(x * 2.3025850929940459f); }
inline double fastpow10(double x) { return exp(x * 2.3025850929940459); }

inline Real logistic(Real x) { return 1.0 / (1.0 + exp(-x)); }
inline Real logit(Real x) { return log(x / (1.0 - x)); }
inline Real dlogit(Real x) { return 1.0 / (x - x*x); }
inline Real logit_scale(Real x, Real a, Real b) { return log((a-x)/(x-b)); }
inline Real dlogit_scale(Real x, Real a, Real b) { return (b-a)/((a-x)*(x-b)); }

inline float rsqrt(float x)
{
	// Based on http://www.hlnum.org/english/doc/frsqrt/frsqrt.html
	__m128 f = _mm_set_ss(x);
	f = _mm_rsqrt_ss(f);
	float r = _mm_cvtss_f32(f);
	r *= ((3.0f - r * r * x) * 0.5f);
	return r;
}

inline double rsqrt(double x)
{
	// Based on http://www.hlnum.org/english/doc/frsqrt/frsqrt.html
	__m128 f = _mm_set_ss((float)x);
	f = _mm_rsqrt_ss(f);
	double r = (double)_mm_cvtss_f32(f);
	r *= ((3.0 - r * r * x) * 0.5);
	r *= ((3.0 - r * r * x) * 0.5);
	return r;
}

inline float logsumf(float loga, float logb)
{
	if (logb > loga) {
		std::swap(loga, logb);
	}
	if (loga == -std::numeric_limits<float>::infinity()) {
		return loga;
	}

	float diff = logb - loga;
	if (diff < -60) {
		return loga;
	}
	float expdiff = exp(diff);
	return loga + boost::math::log1p(expdiff);
}

inline double logsum(double loga, double logb)
{
	if (logb > loga) {
		std::swap(loga, logb);
	}
	if (loga == -std::numeric_limits<double>::infinity()) {
		return loga;
	}

	double diff = logb - loga;
	if (diff < -500) {
		return loga;
	}
	double expdiff = exp(diff);
	return loga + boost::math::log1p(expdiff);
}

inline Real logsum(VectorReal v)
{
	Real m = v.maxCoeff();
	Real sum = 0.0;
	for (unsigned int i = 0; i < v.size(); i++) {
		sum += exp(v(i) - m);
	}
	return log(sum) + m;
}

inline double logaminusb(double loga, double logb)
{
	if (loga == -std::numeric_limits<double>::infinity()) {
		return loga;
	}
	if (logb == -std::numeric_limits<double>::infinity()) {
		return loga;
	}
	if (logb > loga) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	double diff = logb - loga;
	if (diff > 500) {
		return 500;
	} else {
		double expdiff = exp(diff);
		return loga + boost::math::log1p(-expdiff);
	}
}

inline double logabsaminusb(double loga, double logb)
{
	if (logb > loga) {
		std::swap(loga, logb);
	}
	if (loga == -std::numeric_limits<double>::infinity()) {
		return logb;
	}
	if (logb == -std::numeric_limits<double>::infinity()) {
		return loga;
	}

	double diff = logb - loga;
	if (diff > 500) {
		return logb;
	} else {
		double expdiff = exp(diff);
		return loga + boost::math::log1p(-expdiff);
	}
}

inline int factorial(int a)
{
	int result = 1;
	for (int i = 2; i <= a; i++) {
		result *= i;
	}
	return result;
}

inline Real logfactorial_approx(Real a)
{
	//return a * log(a) - a;
	return a * log(a) - a + log(a * (1 + 4 * a * (1 + 2 * a))) / 6.0 + 0.57236494292470008707171367567653;
}

// Copied from https://gist.github.com/jrade
inline float approx_exp(float x)
{
	constexpr float a = (1 << 23) / 0.69314718f;
	constexpr float b = (1 << 23) * (127 - 0.043677448f);
	x = a * x + b;

	constexpr float c = (1 << 23);
	constexpr float d = (1 << 23) * 255;
	if (x < c || x > d)
		x = (x < c) ? 0.0f : d;

	uint32_t n = static_cast<uint32_t>(x);
	memcpy(&x, &n, 4);
	return x;
}

// Copied from https://gist.github.com/jrade
inline double approx_exp(double x)
{
	constexpr double a = (1ll << 52) / 0.6931471805599453;
	constexpr double b = (1ll << 52) * (1023 - 0.04367744890362246);
	x = a * x + b;

	constexpr double c = (1ll << 52);
	constexpr double d = (1ll << 52) * 2047;
	if (x < c || x > d)
		x = (x < c) ? 0.0 : d;

	uint64_t n = static_cast<uint64_t>(x);
	memcpy(&x, &n, 8);
	return x;
}

}
