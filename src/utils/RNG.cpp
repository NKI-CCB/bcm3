#include "Utils.h"
#include "RNG.h"

#if 1

namespace bcm3 {

RNG::RNG(unsigned long long seed)
{
	rl.seed(seed);
}

RNG::~RNG()
{
}

void RNG::Seed(unsigned long long seed)
{
	rl.seed(seed);
}

unsigned int RNG::GetUnsignedInt()
{
	unsigned long long value = rl();
	return (unsigned int)(value & 0xFFFFFFFF);
}

unsigned int RNG::GetUnsignedInt(unsigned int max)
{
	std::uniform_int_distribution<unsigned int> distribution(0, max);
	return distribution(rl);
}

unsigned int RNG::Sample(const VectorReal& probabilities)
{
	if (probabilities.size() == 0) {
		return -1;
	}

	Real t = GetReal();
	Real p = 0.0;
	for (size_t i = 0; i < probabilities.size(); i++) {
		p += probabilities(i);
		if (t < p) {
			return (unsigned int)i;
		}
	}
	return probabilities.size() - 1;
}

Real RNG::GetReal()
{
	return std::generate_canonical<double, std::numeric_limits<double>::digits>(rl);
}

void RNG::fillbuffer() {
#if 0
	uint64_t 
	const uint64_t
		one = 0x3ff0000000000000, // exponent
		m = 0x000fffffffffffff; // mantissa
	uint64_t* id = (uint64_t*)buffer;
	id[0] = one | (m & _x[0]);
	id[1] = one | (m & ((_x[0] >> 52) | (_x[1] << 12)));
	id[2] = one | (m & ((_x[1] >> 40) | (_x[2] << 24)));
	id[3] = one | (m & ((_x[2] >> 28) | (_x[3] << 36)));
	id[4] = one | (m & ((_x[3] >> 16) | (_x[4] << 48)));
	id[5] = one | (m & ((_x[4] >> 4) | (_x[5] << 60)));
	id[6] = one | (m & ((_x[4] >> 56) | (_x[5] << 8)));
	id[7] = one | (m & ((_x[5] >> 44) | (_x[6] << 20)));
	id[8] = one | (m & ((_x[6] >> 32) | (_x[7] << 32)));
	id[9] = one | (m & ((_x[7] >> 20) | (_x[8] << 44)));
	id[10] = one | (m & _x[8] >> 8);
#endif
}

Real RNG::GetGamma(const Real k, const Real theta)
{
	ASSERT(k > 0);

	if (k < 1) {
		const Real u = GetReal();
		return GetGamma((Real)1.0 + k, theta) * pow(u, (Real)1.0 / k);
	} else {
		const Real d = k - 0.33333333333333333333333333333333;
		const Real c = 0.33333333333333333333333333333333 / (sqrt(d));
		Real x, v, u;
		while (1) {
			do {
				x = GetNormal(0.0, 1.0);
				v = 1.0 + c * x;
			} while (v <= 0.0);

			v = v * v * v;
			u = GetReal();

			if (u < 1 - 0.0331 * x * x * x * x)
				break;
			if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
				break;
		}
		return theta * d * v;
	}
}

Real RNG::GetNormal(const Real mu, const Real sigma)
{
	// A bit wasteful
	Real u, v, s;
	do {
		u = GetReal() * 2.0 - 1.0;
		v = GetReal() * 2.0 - 1.0;
		s = u * u + v * v;
	} while (s >= 1.0 || s == 0.0);
	s = sqrt(-2.0 * log(s) / s);
	return mu + sigma * u * s;
}

VectorReal RNG::GetMultivariateUnitNormal(size_t D)
{
	VectorReal res = VectorReal(D);
	for (size_t i = 0; i < D; i++) {
		if (i + 1 < D) {
			Real u, v, s;
			do {
				u = GetReal() * 2.0 - 1.0;
				v = GetReal() * 2.0 - 1.0;
				s = u * u + v * v;
			} while (s >= 1.0 || s == 0.0);
			s = sqrt(-2.0 * log(s) / s);
			res(i) = u * s;
			res(i + 1) = v * s;
			i++;
		} else {
			Real u, v, s;
			do {
				u = GetReal() * 2.0 - 1.0;
				v = GetReal() * 2.0 - 1.0;
				s = u * u + v * v;
			} while (s >= 1.0 || s == 0.0);
			s = sqrt(-2.0 * log(s) / s);
			res(i) = u * s;
		}
	}
	return res;
}

unsigned long long RNG::GetTimeBasedSeed()
{
#if PLATFORM_WINDOWS
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	return (unsigned long long)t.QuadPart;
#elif PLATFORM_LINUX
	timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (unsigned long long)t.tv_nsec;
#elif PLATFORM_MACOSX
    clock_serv_t clock;
    mach_timespec_t t;
    host_get_clock_service(mach_host_self(), REALTIME_CLOCK, &clock);
    clock_get_time(clock, &t);
    mach_port_deallocate(mach_task_self(), clock);
    return (unsigned long long)t.tv_nsec;
#else
	#error
#endif
}
}

#else
#include "../../dependencies/ranluxpp-master/inc/ranluxpp.h"

namespace bcm3 {

RNG::RNG(unsigned long long seed)
{
	rl = new ranluxpp(seed);
}

RNG::~RNG()
{
	delete rl;
}

void RNG::Seed(unsigned long long seed)
{
	delete rl;
	rl = new ranluxpp(seed);
}

unsigned int RNG::GetUnsignedInt()
{
	return rl->getuint();
}

unsigned int RNG::GetUnsignedInt(unsigned int max)
{
	return rl->getuint() % max;
}

unsigned int RNG::Sample(const VectorReal& probabilities)
{
	if (probabilities.size() == 0) {
		return -1;
	}

	Real t = GetReal();
	Real p = 0.0;
	for (size_t i = 0; i < probabilities.size(); i++) {
		p += probabilities(i);
		if (t < p) {
			return (unsigned int)i;
		}
	}
	return probabilities.size() - 1;
}

Real RNG::GetReal()
{
	return rl->getdouble();
}

Real RNG::GetGamma(const Real k, const Real theta)
{
	ASSERT(k > 0);

	if (k < 1) {
		const Real u = GetReal();
		return GetGamma((Real)1.0 + k, theta) * pow(u, (Real)1.0 / k);
	} else {
		const Real d = k - 0.33333333333333333333333333333333;
		const Real c = 0.33333333333333333333333333333333 / (sqrt(d));
		Real x, v, u;
		while (1) {
			do {
				x = GetNormal(0.0, 1.0);
				v = 1.0 + c * x;
			} while (v <= 0.0);

			v = v * v * v;
			u = GetReal();

			if (u < 1 - 0.0331 * x * x * x * x)
				break;
			if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
				break;
		}
		return theta * d * v;
	}
}

Real RNG::GetNormal(const Real mu, const Real sigma)
{
	// A bit wasteful
	Real u, v, s;
	do {
		u = GetReal() * 2.0 - 1.0;
		v = GetReal() * 2.0 - 1.0;
		s = u * u + v * v;
	} while (s >= 1.0 || s == 0.0);
	s = sqrt(-2.0 * log(s) / s);
	return mu + sigma * u * s;
}

unsigned long long RNG::GetTimeBasedSeed()
{
#if PLATFORM_WINDOWS
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	return (unsigned long long)t.QuadPart;
#elif PLATFORM_LINUX
	timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (unsigned long long)t.tv_nsec;
#elif PLATFORM_MACOSX
    clock_serv_t clock;
    mach_timespec_t t;
    host_get_clock_service(mach_host_self(), REALTIME_CLOCK, &clock);
    clock_get_time(clock, &t);
    mach_port_deallocate(mach_task_self(), clock);
    return (unsigned long long)t.tv_nsec;
#else
	#error
#endif
}
}
#endif

