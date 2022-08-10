#include "Utils.h"
#include "Logger.h"

#include <stdarg.h>
#include <boost/date_time.hpp>
using namespace boost::posix_time;

namespace bcm3 {

Logger* logger = NULL;

class LoggerInitializer
{
public:
	LoggerInitializer()
	{
		logger = new Logger;
	}
	~LoggerInitializer()
	{
		delete logger;
	}
};

LoggerInitializer init;

static const int MAX_MESSAGE_SIZE = 16 * 1024;
static const char LEVEL_TO_STRING_MAP[Logger::Count] = { ' ', 'W', 'E', '?' };

Logger::Logger()
: File(NULL)
, MessageBuffer(new char[MAX_MESSAGE_SIZE])
, ConsoleLevel(None)
, FileLevel(None)
, BufferLevel(None)
, LogStartTime(second_clock::local_time())
{
}

Logger::~Logger()
{
	if (File) {
		fclose(File);
	}
	delete[] MessageBuffer;
}

bool Logger::SetLogToConsole(ELevel level)
{
	std::lock_guard<std::mutex> lock(LogMutex);

	ConsoleLevel = level;
	return true;
}

bool Logger::SetLogToFile(ELevel level, const char* filename)
{
	std::lock_guard<std::mutex> lock(LogMutex);

	if (File) {
		fclose(File);
		File = NULL;
	}
	FileLevel = None;

	if (level != None) {
		File = fopen(filename, "w");
		if (File) {
			FileLevel = level;
			return true;
		} else {
			return false;
		}
	} else {
		return true;
	}
}

void Logger::Log(ELevel level, const char* source, int line, const char* message, ...)
{
	std::lock_guard<std::mutex> lock(LogMutex);

	if (level < Info || level >= Count) {
		return;
	}

	unsigned long threadid;
#if PLATFORM_WINDOWS
	threadid = GetCurrentThreadId();
	const char* short_source = strrchr(source, '\\') ? strrchr(source, '\\') + 1 : source;
#elif PLATFORM_LINUX
	threadid = pthread_self();
	const char* short_source = strrchr(source, '/') ? strrchr(source, '/') + 1 : source;
#elif PLATFORM_MACOSX
	threadid = (unsigned long)pthread_mach_thread_np(pthread_self());
	const char* short_source = strrchr(source, '/') ? strrchr(source, '/') + 1 : source;
#else
	#error
#endif
	
	time_duration time = second_clock::local_time() - LogStartTime;
	snprintf(MessageBuffer, MAX_MESSAGE_SIZE, "%1.1s %04lx %02u:%02u:%02u %-16.16s:%-5d| ",
		&LEVEL_TO_STRING_MAP[level],
		(threadid % 0xffff),
		(unsigned int)time.hours(),
		(unsigned int)time.minutes(),
		(unsigned int)time.seconds(),
		short_source,
		line);
	
	va_list args;
	va_start(args, message);
	vsnprintf(MessageBuffer + 40, MAX_MESSAGE_SIZE - 40, message, args);
	va_end(args);

	LogConsole(level, MessageBuffer);
	LogFile(level, MessageBuffer);
}

void Logger::LogTime()
{
	boost::posix_time::ptime time = boost::date_time::second_clock<boost::posix_time::ptime>::local_time();
	LOG("Current date: %04u-%02u-%02u",
			(unsigned int)time.date().year(),
			(unsigned int)time.date().month(),
			(unsigned int)time.date().day());
	LOG("Current time: %02u:%02u:%02u",
			(unsigned int)time.time_of_day().hours(),
			(unsigned int)time.time_of_day().minutes(),
			(unsigned int)time.time_of_day().seconds());
}

void Logger::LogCPUInfo()
{
#if PLATFORM_LINUX
	FILE* file = fopen("/proc/cpuinfo", "r");
	if (file) {
		char buf[1024];
		while(fgets(buf, 1024, file)) {
			size_t newline = strlen(buf)-1;
			buf[newline] = 0;
			if (strstr(buf, "model name") != NULL) {
				LOG(buf);
			}
		}
		fclose(file);
	}
#else
	LOG("(Only implemented on Linux)");
#endif
}

void Logger::LogConsole(ELevel level, const char* message)
{
	if (level >= ConsoleLevel) {
		puts(message);
	}
}

void Logger::LogFile(ELevel level, const char* message)
{
	if (level >= FileLevel) {
		fputs(message, File);
		fputs("\n", File);
		fflush(File);
	}
}

}
