#include "Utils.h"
#include "ProgressIndicatorConsole.h"

#if PLATFORM_LINUX
#include <unistd.h>
#endif

namespace bcm3 {

ProgressIndicatorConsole::ProgressIndicatorConsole()
	: prevlen(0)
	, previous_time(-std::numeric_limits<Real>::infinity())
	, update_settings(0.5)
{
}

ProgressIndicatorConsole::~ProgressIndicatorConsole()
{
}

void ProgressIndicatorConsole::SetUpdateTime(double seconds)
{
	update_settings = seconds;
}

void ProgressIndicatorConsole::OnUpdate(bool force)
{
	// Don't update more than once every half second
	if (previous_time == -std::numeric_limits<Real>::infinity()) {
		timer.Start();
		previous_time = timer.GetElapsedSeconds();
	} else {
		double new_time = timer.GetElapsedSeconds();
		if (new_time - previous_time > update_settings || current_progress >= 1.0) {
			previous_time = new_time;
		} else if (!force) {
			return;
		}
	}

	char buf[256];
	for (size_t i = 0; i < prevlen; i++) {
		buf[i] = ' ';
	}
	buf[prevlen] = 0;

#if PLATFORM_LINUX
	if (isatty(fileno(stdout))) {
		printf("\r%s", buf);
	}
#else
	printf("\r%s", buf);
#endif

	snprintf(buf, sizeof(buf), "Progress: %6.2f%% - remaining: %s", current_progress*100.0, eta.c_str());
	for (std::map<std::string, std::string>::iterator info = current_info.begin(); info != current_info.end(); ++info) {
		strncat(buf, " - ", sizeof(buf) - strlen(buf) - 1);
		strncat(buf, info->first.c_str(), sizeof(buf) - strlen(buf) - 1);
		strncat(buf, ": ", sizeof(buf) - strlen(buf) - 1);
		strncat(buf, info->second.c_str(), sizeof(buf) - strlen(buf) - 1);
	}

#if PLATFORM_LINUX
	if (isatty(fileno(stdout))) {
		printf("\r%s", buf);
	} else {
		printf("%s\n", buf);
	}
#else
	printf("\r%s", buf);
#endif
	fflush(stdout);

	prevlen = strlen(buf);
}

void ProgressIndicatorConsole::OnStop()
{
	printf("\n");
}

}
