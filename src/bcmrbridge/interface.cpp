#include "Utils.h"
#include "interface.h"

#include <boost/filesystem.hpp>

bcm3info* GetBCM3InfoPtr(char** bcm3info_ptr, int* retval)
{
	void* p;
	int res = sscanf(*bcm3info_ptr, "%p", &p);
	if (res != 1) {
		LOGERROR("Could not retrieve bcm3 info pointer");
		*retval = -1;
		return NULL;
	}
	bcm3info* info = (bcm3info*)p;
	if (!info) {
		LOGERROR("Invalid bcm3 info pointer");
		*retval = -1;
		return NULL;
	} else {
		return info;
	}
}

extern "C" {

void bcm3_rbridge_init(char** bcm3info_ptr, char** base_folder, char** prior_fn, char** likelihood_fn, int* num_threads, int* retval)
{
	bcm3info* info = new bcm3info;
	if (!info) {
		*retval = -1;
		return;
	}

	bcm3::logger->SetLogToFile(bcm3::Logger::Info, "bcm3rbridge.log");

	boost::filesystem::path cwd = boost::filesystem::current_path();
	boost::filesystem::current_path(*base_folder);

	info->varset = std::make_shared<bcm3::VariableSet>();
	if (!info->varset->LoadFromXML(*prior_fn)) {
		*retval = -2;
		boost::filesystem::current_path(cwd);
		return;
	}

	info->prior = bcm3::Prior::Create(*prior_fn, info->varset, 1);
	if (!info->prior) {
		*retval = -3;
		boost::filesystem::current_path(cwd);
		return;
	}

	size_t evaluation_threads = 1;
	if (*num_threads == -1) {
		evaluation_threads = std::max(1u, std::thread::hardware_concurrency() - 1);
	} else if (*num_threads > 0) {
		evaluation_threads = (size_t)*num_threads;
	} else {
		LOGERROR("Invalid number of threads");
		*retval = -4;
		return;
	}

	info->likelihood = LikelihoodFactory::CreateLikelihood(*likelihood_fn, info->varset, 1, evaluation_threads);
	if (!info->likelihood) {
		*retval = -5;
		boost::filesystem::current_path(cwd);
		return;
	}
	info->likelihood->PostInitialize();
	
	snprintf(*bcm3info_ptr, 128, "%p", (void*)info);
	boost::filesystem::current_path(cwd);
	*retval = 0;
}

}
