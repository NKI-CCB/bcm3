#include "Utils.h"
#include "BlockingStrategy.h"

namespace bcm3 {

	void BlockingStrategy::Initialize(size_t num_variables)
	{
		this->num_variables = num_variables;
	}

}
