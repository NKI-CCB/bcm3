#pragma once

#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "Version.h"
#include "Platform.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#if PLATFORM_MACOSX
	#include <malloc/malloc.h>
#else
	#include <malloc.h>
#endif

#include <vector>
#include <list>
#include <set>
#include <queue>
#include <map>
#include <algorithm>
#include <limits>
#include <string>
#include <functional>
#include <thread>
#include <mutex>
#include <memory>

#include <boost/version.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>

#include "Debug.h"

#if PLATFORM_WINDOWS && BUILD_DEBUG
	#define EIGEN_RUNTIME_NO_MALLOC
#endif
#include "Eigen/Dense"

#include "Typedefs.h"
#include "Logger.h"
#include "MathFunctions.h"

namespace bcm3 {

void tokenize(std::string str, std::vector<std::string>& tokens, const char* delims);
void fix_path(std::string& path);

}
