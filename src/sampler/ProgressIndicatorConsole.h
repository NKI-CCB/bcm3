#pragma once

#include "ProgressIndicator.h"

namespace bcm3 {

class ProgressIndicatorConsole : public ProgressIndicator
{
public:
	ProgressIndicatorConsole();
	virtual ~ProgressIndicatorConsole();

	void SetUpdateTime(double seconds);

private:
	virtual void OnUpdate(bool force);
	virtual void OnStop();

	size_t prevlen;
	Timer timer;
	double previous_time;
	double update_settings;
};

}
