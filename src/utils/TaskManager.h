#pragma once

#include <condition_variable>

namespace bcm3 {

typedef std::function<bool (void*, size_t)> TTask;

class TaskManager
{
public:
	TaskManager(size_t numthreads);
	~TaskManager();

	uint64 AddTask(TTask task, void* userdata);
	bool WaitTask(uint64 task);

	inline size_t GetNumThreads() const { return Threads.size(); }

private:
	struct Task
	{
		TTask Delegate;
		void* UserData;
		uint64 Handle;
		std::atomic<bool> return_value;
	};
	void WorkerFunction(size_t threadIndex);
	friend void worker(TaskManager* manager, size_t threadIndex);

	std::mutex StartMutex;
	std::mutex FinishMutex;
	std::condition_variable StartCondition;
	std::condition_variable FinishCondition;
	bool Active;

	std::vector<std::thread*> Threads;
	std::list<Task*> WaitingTasks;
	std::list<Task*> FinishedTasks;
	uint64 NextHandle;
};

}
