module rm @private @feat(RM_THREADED);

import std::atomic::types;
import std::thread;
import std::math;
import std::os;

// ================================
// Private: Threads Functions
// ================================

fn void threads_init(int total)
{
	if (total <= 0) total = os::num_cpu();

	ctx.workers = math::max(total - 1, 0);
	ctx.job.mx.init()!!;
	ctx.job.done.init()!!;

	if (ctx.workers > 0)
	{
		ctx.pool.init(mem, ctx.workers)!!;
		ctx.pool.prewarm(ctx.workers)!!;
	}
}

fn void threads_quit()
{
	if (ctx.workers > 0) ctx.pool.destroy();

	ctx.job.done.destroy();
	ctx.job.mx.destroy();
}

fn void threads_run(ItemFn item_fn, void* arg, int count)
{
	Job* job = &ctx.job;

	job.item_fn = item_fn;
	job.arg     = arg;
	job.count   = count;
	job.next.store(0, RELAXED);

	int helpers = math::max(math::min(ctx.workers, count - 1), 0);

	job.mx.lock();
	job.pending = helpers;
	job.mx.unlock();

	for (int i = 0; i < helpers; i++)
	{
		if (catch ctx.pool.add_task(&worker_main, job))
		{
			job_finish(job);
		}
	}

	run_items(job);

	job.mx.lock();
	while (job.pending > 0) job.done.wait(&job.mx);
	job.mx.unlock();
}

// ================================
// Local: Threads Functions
// ================================

fn void worker_main(void* arg) @local
{
	Job* job = arg;

	run_items(job);
	job_finish(job);
}

// Claims and runs items until there are none left
fn void run_items(Job* job) @local
{
	while (true)
	{
		int item = job.next.add(1, RELAXED); // returns the previous value
		if (item >= job.count) break;
		job.item_fn(job.arg, item);
	}
}

fn void job_finish(Job* job) @local
{
	job.mx.lock();
	job.pending--;
	if (job.pending == 0) job.done.signal();
	job.mx.unlock();
}

// ================================
// Local: Threads Context
// ================================

struct Job @local
{
	// Current run (read-only)
	ItemFn item_fn;
	void*  arg;
	int    count;

	Atomic{int} next; // next item to claim

	// Completion
	Mutex             mx;
	ConditionVariable done;
	int               pending; // workers that haven't finished (protected by mx)
}

struct Context @local
{
	ThreadPool pool;
	Job        job;
	int        workers;
}

Context ctx @local;
