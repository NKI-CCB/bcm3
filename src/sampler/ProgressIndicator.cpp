#include "Utils.h"
#include "ProgressIndicator.h"

namespace bcm3 {

ProgressIndicator::ProgressIndicator()
	: current_run_time(0.0)
	, expected_run_time(0.0)
	, expected_remaining_time(0.0)
	, expected_remaining_time_ema(0.0)
	, min_update_time(std::numeric_limits<Real>::quiet_NaN())
	, min_update_progress(std::numeric_limits<Real>::quiet_NaN())
	, last_progress(-std::numeric_limits<Real>::max())
	, last_time(-std::numeric_limits<Real>::max())
{
}

ProgressIndicator::~ProgressIndicator()
{
}

void ProgressIndicator::NotifyStart()
{
	timer.Start();
}

void ProgressIndicator::NotifyStop()
{
	OnStop();
}

void ProgressIndicator::UpdateInfo(const char* name, const std::string& value)
{
	current_info[name] = value;
}

void ProgressIndicator::SetMinUpdateTime(Real t)
{
	min_update_time = t;
}

void ProgressIndicator::SetMinUpdateProgress(Real progress)
{
	min_update_progress = progress;
}

void ProgressIndicator::ReportProgress(Real progress, bool force)
{
	current_progress = progress;
	current_run_time = std::max(1e-5, timer.GetElapsedSeconds());
	if (progress <= 0) {
		progress = 1e-5;
	}
	expected_run_time = current_run_time / progress;
	expected_remaining_time = current_run_time * (1.0 - progress) / progress;
	expected_remaining_time_ema += (expected_remaining_time - expected_remaining_time_ema) * 0.1;

	if (expected_remaining_time > (Real)std::numeric_limits<uint64>::max()) {
		eta = "infinite";
	} else {
		uint64 seconds = (uint64)expected_remaining_time;
		uint64 minutes = seconds / 60;
		seconds -= minutes * 60;
		uint64 hours = minutes / 60;
		minutes -= hours * 60;
		uint64 days = hours / 24;
		hours -= days * 24;

		char buf[32];
		snprintf(buf, sizeof(buf), "%llud %02lluh %02llum %02llus", days, hours, minutes, seconds);
		eta = buf;
	}

	bool do_update = true;
	if (min_update_time == min_update_time) {
		if (current_run_time < last_time + min_update_time) {
			do_update = false;
		}
	}
	if (min_update_progress == min_update_progress) {
		if (progress < last_progress + min_update_progress) {
			do_update = false;
		}
	}
	if (progress >= 1.0) {
		do_update = true;
	}
	if (do_update || force) {
		OnUpdate(force);
		last_progress = progress;
		last_time = current_run_time;
	}
}

}
