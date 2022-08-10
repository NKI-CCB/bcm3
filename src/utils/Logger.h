#pragma once

#include <boost/date_time/posix_time/posix_time.hpp>

#define LOG(message, ...)			do { bcm3::logger->Log(bcm3::Logger::Info   , __FILE__, __LINE__, message, ## __VA_ARGS__); } while (0)
#define LOGWARNING(message, ...)	do { bcm3::logger->Log(bcm3::Logger::Warning, __FILE__, __LINE__, message, ## __VA_ARGS__); } while (0)
#define LOGERROR(message, ...)		do { bcm3::logger->Log(bcm3::Logger::Error  , __FILE__, __LINE__, message, ## __VA_ARGS__); } while (0)

namespace bcm3 {

class Logger
{
public:
	enum ELevel
	{
		Info	= 0,
		Warning = 1,
		Error	= 2,
		None	= 3,

		Count
	};

	Logger();
	~Logger();

	bool SetLogToConsole(ELevel level);
	bool SetLogToFile(ELevel level, const char* filename);

	void Log(ELevel level, const char* source_file, int line, const char* message, ...);
	void LogTime();
	void LogCPUInfo();

private:
	void LogConsole(ELevel level, const char* message);
	void LogFile(ELevel level, const char* message);

	FILE* File;
	char* MessageBuffer;

	ELevel ConsoleLevel;
	ELevel FileLevel;
	ELevel BufferLevel;

	boost::posix_time::ptime LogStartTime;

	std::mutex LogMutex;
};

extern Logger* logger;

}
