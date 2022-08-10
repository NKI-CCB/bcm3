#pragma once

#ifdef _WIN32
	#define PLATFORM_WINDOWS 1
	#define NOMINMAX
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
#else
	#define PLATFORM_WINDOWS 0
#endif

#ifdef __linux__
	#define PLATFORM_LINUX 1
#else
	#define PLATFORM_LINUX 0
#endif

#if defined(__APPLE__) && defined(__MACH__)
	// Could be iOS as well I think, but close enough
	#define PLATFORM_MACOSX 1
#else
	#define PLATFORM_MACOSX 0
#endif

#if PLATFORM_WINDOWS
	#if _MSC_VER < 1700
		#define snprintf _snprintf
	#endif
	#define strcasecmp _stricmp
#elif PLATFORM_LINUX
	#define _aligned_malloc(size, alignment) memalign(alignment, size)
	#define _aligned_free(p) free(p)
#elif PLATFORM_MACOSX
	#define _aligned_malloc(size, alignment) malloc(size)
	#define _aligned_free(p) free(p)
#endif
