#pragma once

#include "SampleHandler.h"

namespace bcm3 {

	class SampleHandlerStoreMaxAPosteriori: public SampleHandler
	{
	public:
		SampleHandlerStoreMaxAPosteriori();
		~SampleHandlerStoreMaxAPosteriori();

		void Reset();
	
		virtual void ReceiveSample(const VectorReal& values, Real lprior, Real llh, Real temperature, Real weight);

		inline const Real GetMAPlposterior() const { return MAP_lposterior; }
		inline const VectorReal& GetMAP() const { return MAP_values; }

	private:
		Real MAP_lposterior;
		VectorReal MAP_values;
	};

}
