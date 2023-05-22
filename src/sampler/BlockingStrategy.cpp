#include "Utils.h"
#include "BlockingStrategy.h"

namespace bcm3 {

	BlockingStrategy::BlockingStrategy()
		: num_variables(0)
	{
	}

	BlockingStrategy::~BlockingStrategy()
	{
	}

	void BlockingStrategy::Initialize(size_t num_variables)
	{
		this->num_variables = num_variables;
	}

	bool BlockingStrategy::UsesClustering()
	{
		return false;
	}
}
