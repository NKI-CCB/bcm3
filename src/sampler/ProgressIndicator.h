#pragma once

#include "Timer.h"

namespace bcm3 {

class ProgressIndicator
{
public:
	ProgressIndicator();
	virtual ~ProgressIndicator();

	void NotifyStart();
	void NotifyStop();
	void UpdateInfo(const char* name, const std::string& value);

	void SetMinUpdateTime(Real t);
	void SetMinUpdateProgress(Real progress);

	/// Call UpdateInfo before ReportProgress
	void ReportProgress(Real progress, bool force = false);

protected:
	virtual void OnUpdate(bool force) {}
	virtual void OnStop() {}

	Real current_progress;
	std::map<std::string, std::string> current_info;

	Timer timer;
	Real current_run_time;
	Real expected_run_time;
	Real expected_remaining_time;
	Real expected_remaining_time_ema;
	std::string eta;

	Real min_update_time;
	Real min_update_progress;
	Real last_progress;
	Real last_time;
};

}
