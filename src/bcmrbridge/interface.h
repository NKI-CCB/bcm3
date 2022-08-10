#pragma once

#include "LikelihoodFactory.h"
#include "Prior.h"
#include "VariableSet.h"

struct bcm3info {
	std::shared_ptr<bcm3::Prior> prior;
	std::shared_ptr<bcm3::Likelihood> likelihood;
	std::shared_ptr<bcm3::VariableSet> varset;
};
bcm3info* GetBCM3InfoPtr(char** bcm3info_ptr, int* retval);
