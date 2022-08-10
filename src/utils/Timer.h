#pragma once

namespace bcm3 {

#if PLATFORM_WINDOWS

class Timer
{
public:
	void Start()
	{
		QueryPerformanceFrequency(&Frequency);
		QueryPerformanceCounter(&StartTime);
	}

	double GetElapsedSeconds()
	{
		QueryPerformanceCounter(&EndTime);
		double elapsed = (EndTime.QuadPart - StartTime.QuadPart) / (double)(Frequency.QuadPart);
		return elapsed;
	}

private:
	LARGE_INTEGER Frequency;
	LARGE_INTEGER StartTime, EndTime;
};

#elif PLATFORM_LINUX

#include <time.h>

class Timer
{
public:
	void Start()
	{
		clock_gettime(CLOCK_MONOTONIC, &StartTime);
	}

	double GetElapsedSeconds()
	{
		timespec end;
		clock_gettime(CLOCK_MONOTONIC, &end);
		double elapsed = (end.tv_sec - StartTime.tv_sec) + (end.tv_nsec - StartTime.tv_nsec) * 0.000000001;
		return elapsed;
	}

private:
	timespec StartTime;
};

#elif PLATFORM_MACOSX
	
#include <mach/clock.h>
#include <mach/mach.h>
    
    class Timer
    {
    public:
		Timer()
		{
			host_get_clock_service(mach_host_self(), REALTIME_CLOCK, &Clock);
		}

		~Timer()
		{
			mach_port_deallocate(mach_task_self(), Clock);
		}

        void Start()
        {
			clock_get_time(Clock, &StartTime);
        }
        
        double GetElapsedSeconds()
        {
			mach_timespec_t end;
			clock_get_time(Clock, &end);
			double elapsed = (end.tv_sec - StartTime.tv_sec) + (end.tv_nsec - StartTime.tv_nsec) * 0.000000001;
			return elapsed;
        }
        
    private:
		mach_timespec_t StartTime;
		clock_serv_t Clock;
    };

    
#endif

}
