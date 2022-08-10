#include "Utils.h"
#include "TaskManager.h"
#include <boost/interprocess/detail/atomic.hpp>

namespace bcm3 {

void worker(TaskManager* manager, size_t threadIndex)
{
	manager->WorkerFunction(threadIndex);
}

TaskManager::TaskManager(size_t numthreads)
	: NextHandle(0)
	, Active(true)
{
	Threads.resize(numthreads);
	for (size_t i = 0; i < numthreads; i++) {
		Threads[i] = new std::thread(worker, this, i);
	}
}

TaskManager::~TaskManager()
{
	{
		std::lock_guard<std::mutex> lock(StartMutex);
		Active = false;
		StartCondition.notify_all();
	}
	
	for (size_t i = 0; i < Threads.size(); i++) {
		Threads[i]->join();
		delete Threads[i];
	}
}

uint64 TaskManager::AddTask(TTask task, void* userdata)
{
	std::lock_guard<std::mutex> lock(StartMutex);
	Task* t = new Task;
	t->Delegate = task;
	t->UserData = userdata;
	t->Handle = NextHandle;
	t->return_value = false;
	NextHandle++;
	WaitingTasks.push_back(t);
	StartCondition.notify_one();
	return t->Handle;
}

bool TaskManager::WaitTask(uint64 task)
{
	bool return_value = false;

	std::unique_lock<std::mutex> lock(FinishMutex);
    while(1) {
		bool finished = false;
		for (std::list<Task*>::iterator ti = FinishedTasks.begin(); ti != FinishedTasks.end(); ++ti) {
			if ((*ti)->Handle == task) {
				finished = true;
				return_value = (*ti)->return_value.load(std::memory_order_acquire);
				delete *ti;
				FinishedTasks.erase(ti);
				break;
			}
		}

		if (finished) {
			break;
		} else {
			FinishCondition.wait(lock);
		}
    }

	return return_value;
}

void TaskManager::WorkerFunction(size_t threadIndex)
{
	while (1) {
		Task* task;

		// Wait for a task.
		{
			std::unique_lock<std::mutex> lock(StartMutex);
			std::list<Task*>::iterator ti;
			while (1) {
				if (!Active) {
					break;
				}
				ti = WaitingTasks.begin();
				if (ti != WaitingTasks.end()) {
					task = *ti;
					WaitingTasks.erase(ti);
					break;
				} else {
					StartCondition.wait(lock);
				}
			}
		}
			
		if (!Active) {
			break;
		}

		// Execute task.
		bool retval = task->Delegate(task->UserData, threadIndex);

		// Append task to finished lists.
		{
			std::lock_guard<std::mutex> lock(FinishMutex);
			FinishedTasks.push_back(task);
			task->return_value.store(retval, std::memory_order_release);
		}
		FinishCondition.notify_one();
	}
}

}
