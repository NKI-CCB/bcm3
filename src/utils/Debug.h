#pragma once

#ifdef _DEBUG
	#define BUILD_DEBUG 1
	#define BUILD_RELEASE 0
#else
	#define BUILD_DEBUG 0
	#define BUILD_RELEASE 1
#endif

#if BUILD_DEBUG
	#if PLATFORM_WINDOWS
		#define ASSERT(expr)																		\
			do																						\
			{																						\
				if (!(expr))																		\
				{																					\
					bcm3::logger->Log(bcm3::Logger::Error, __FILE__, __LINE__, "assert " , #expr);	\
					__debugbreak();																	\
				}																					\
			} while (0)
	#else
		#error
	#endif
#else
	#define ASSERT(expr) do {} while (0)
#endif
