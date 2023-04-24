#pragma once

namespace bcm3 {

	class SampleHandler
	{
	public:
		virtual void ReceiveSample(const VectorReal& sample, Real lprior, Real llh, Real temperature, Real weight) = 0;
	};

}
