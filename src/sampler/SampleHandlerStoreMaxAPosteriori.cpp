#include "Utils.h"
#include "SampleHandlerStoreMaxAPosteriori.h"

namespace bcm3 {

SampleHandlerStoreMaxAPosteriori::SampleHandlerStoreMaxAPosteriori()
	: MAP_lposterior(-std::numeric_limits<Real>::infinity())
{
}

void SampleHandlerStoreMaxAPosteriori::Reset()
{
	MAP_lposterior = -std::numeric_limits<Real>::infinity();
}

SampleHandlerStoreMaxAPosteriori::~SampleHandlerStoreMaxAPosteriori()
{
}

void SampleHandlerStoreMaxAPosteriori::ReceiveSample(const VectorReal& values, Real lprior, Real llh, Real temperature)
{
	// Use sample regardless of temperature
	Real lposterior = lprior + llh;
	if (lposterior > MAP_lposterior) {
		MAP_lposterior = lposterior;
		MAP_values = values;
	}
}
}
