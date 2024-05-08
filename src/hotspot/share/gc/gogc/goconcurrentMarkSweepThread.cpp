#include "precompiled.hpp"
#include "classfile/systemDictionary.hpp"
#include "gc/gogc/gogcHeap.hpp"
#include "gc/gogc/concurrentMarkSweepGeneration.inline.hpp"
#include "gc/gogc/goConcurrentMarkSweepThread.hpp"
#include "gc/shared/gcId.hpp"
#include "memory/universe.hpp"
#include "oops/oop.inline.hpp"
#include "runtime/init.hpp"
#include "runtime/java.hpp"
#include "runtime/javaCalls.hpp"
#include "runtime/mutexLocker.hpp"
#include "runtime/os.hpp"
#include "runtime/vmThread.hpp"


goConcurrentMarkSweepThread* goConcurrentMarkSweepThread::_gogct = NULL;
gogcCollector* goConcurrentMarkSweepThread::_collector         = NULL;
int  goConcurrentMarkSweepThread::_gogc_flag                   = gogc_nil;

volatile jint goConcurrentMarkSweepThread::_pending_yields    = 0;

goConcurrentMarkSweepThread::goConcurrentMarkSweepThread(gogcCollector* collector)
  : ConcurrentGCThread() {
  assert(_gogct == NULL, "gogc thread already created");
  _gogct = this;
  assert(_collector == NULL, "Collector already set");
  _collector = collector;

  set_name("gogc Main Thread");

  create_and_start(UseCriticalgogcThreadPriority ? CriticalPriority : NearMaxPriority);
}

void goConcurrentMarkSweepThread::run_service() {
  assert(this == gogct(), "just checking");

  if (BindgogcThreadToCPU && !os::bind_to_processor(CPUForgogcThread)) {
    log_warning(gc)("Couldn't bind gogc thread to processor " UINTX_FORMAT, CPUForgogcThread);
  }

  while (!should_terminate()) {
    sleepBeforeNextCycle();
    if (should_terminate()) break;
    GCIdMark gc_id_mark;
    GCCause::Cause cause = _collector->_full_gc_requested ?
      _collector->_full_gc_cause : GCCause::_gogc_concurrent_mark;
    _collector->collect_in_background(cause);
  }

  verify_ok_to_terminate();
}

#ifndef PRODUCT
void goConcurrentMarkSweepThread::verify_ok_to_terminate() const {
  assert(!(CGC_lock->owned_by_self() || gogc_thread_has_gogc_token() ||
           gogc_thread_wants_gogc_token()),
         "Must renounce all worldly possessions and desires for nirvana");
  _collector->verify_ok_to_terminate();
}
#endif

// create and start a new ConcurrentMarkSweep Thread for given gogc generation
goConcurrentMarkSweepThread* goConcurrentMarkSweepThread::start(gogcCollector* collector) {
  guarantee(_gogct == NULL, "start() called twice!");
  goConcurrentMarkSweepThread* th = new goConcurrentMarkSweepThread(collector);
  assert(_gogct == th, "Where did the just-created gogc thread go?");
  return th;
}

void goConcurrentMarkSweepThread::stop_service() {
  MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
  CGC_lock->notify_all();
}

void goConcurrentMarkSweepThread::threads_do(ThreadClosure* tc) {
  assert(tc != NULL, "Null ThreadClosure");
  if (gogct() != NULL && !gogct()->has_terminated()) {
    tc->do_thread(gogct());
  }
  assert(Universe::is_fully_initialized(),
  if (_collector != NULL) {
    AbstractWorkGang* gang = _collector->conc_workers();
    if (gang != NULL) {
      gang->threads_do(tc);
    }
  }
}

void goConcurrentMarkSweepThread::print_all_on(outputStream* st) {
  if (gogct() != NULL && !gogct()->has_terminated()) {
    gogct()->print_on(st);
    st->cr();
  }
  if (_collector != NULL) {
    AbstractWorkGang* gang = _collector->conc_workers();
    if (gang != NULL) {
      gang->print_worker_threads_on(st);
    }
  }
}

void goConcurrentMarkSweepThread::synchronize(bool is_gogc_thread) {
  assert(UseConcMarkSweepGC, "just checking");

  MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
  if (!is_gogc_thread) {
    assert(Thread::current()->is_VM_thread(), "Not a VM thread");
    gogcSynchronousYieldRequest yr;
    while (gogc_flag_is_set(gogc_gogc_has_token)) {
      // indicate that we want to get the token
      set_gogc_flag(gogc_vm_wants_token);
      CGC_lock->wait_without_safepoint_check();
    }
    // claim the token and proceed
    clear_gogc_flag(gogc_vm_wants_token);
    set_gogc_flag(gogc_vm_has_token);
  } else {
    assert(Thread::current()->is_ConcurrentGC_thread(),
           "Not a gogc thread");
    while (gogc_flag_is_set(gogc_vm_has_token | gogc_vm_wants_token)) {
      set_gogc_flag(gogc_gogc_wants_token);
      CGC_lock->wait_without_safepoint_check();
    }
    // claim the token
    clear_gogc_flag(gogc_gogc_wants_token);
    set_gogc_flag(gogc_gogc_has_token);
  }
}

void goConcurrentMarkSweepThread::desynchronize(bool is_gogc_thread) {
  assert(UseConcMarkSweepGC, "just checking");

  MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
  if (!is_gogc_thread) {
    assert(Thread::current()->is_VM_thread(), "Not a VM thread");
    assert(gogc_flag_is_set(gogc_vm_has_token), "just checking");
    clear_gogc_flag(gogc_vm_has_token);
    if (gogc_flag_is_set(gogc_gogc_wants_token)) {
      // wake-up a waiting gogc thread
      CGC_lock->notify();
    }
    assert(!gogc_flag_is_set(gogc_vm_has_token | gogc_vm_wants_token),
           "Should have been cleared");
  } else {
    assert(Thread::current()->is_ConcurrentGC_thread(),
           "Not a gogc thread");
    assert(gogc_flag_is_set(gogc_gogc_has_token), "just checking");
    clear_gogc_flag(gogc_gogc_has_token);
    if (gogc_flag_is_set(gogc_vm_wants_token)) {
      // wake-up a waiting VM thread
      CGC_lock->notify();
    }
    assert(!gogc_flag_is_set(gogc_gogc_has_token | gogc_gogc_wants_token),
           "Should have been cleared");
  }
}

// Wait until any gogc_lock event
void goConcurrentMarkSweepThread::wait_on_gogc_lock(long t_millis) {
  MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
  if (should_terminate() || _collector->_full_gc_requested) {
    return;
  }
  set_gogc_flag(gogc_gogc_wants_token);   // to provoke notifies
  CGC_lock->wait_without_safepoint_check(t_millis);
  clear_gogc_flag(gogc_gogc_wants_token);
  assert(!gogc_flag_is_set(gogc_gogc_has_token | gogc_gogc_wants_token),
         "Should not be set");
}


void goConcurrentMarkSweepThread::wait_on_gogc_lock_for_scavenge(long t_millis) {
  // Wait time in millis or 0 value representing infinite wait for a scavenge
  assert(t_millis >= 0, "Wait time for scavenge should be 0 or positive");

  gogcHeap* heap = gogcHeap::heap();
  double start_time_secs = os::elapsedTime();
  double end_time_secs = start_time_secs + (t_millis / ((double) MILLIUNITS));

  // Total collections count before waiting loop
  unsigned int before_count;
  {
    MutexLocker hl(Heap_lock, Mutex::_no_safepoint_check_flag);
    before_count = heap->total_collections();
  }

  unsigned int loop_count = 0;

  while(!should_terminate()) {
    double now_time = os::elapsedTime();
    long wait_time_millis;

    if(t_millis != 0) {
      // New wait limit
      wait_time_millis = (long) ((end_time_secs - now_time) * MILLIUNITS);
      if(wait_time_millis <= 0) {
        // Wait time is over
        break;
      }
    } else {
      // No wait limit, wait if necessary forever
      wait_time_millis = 0;
    }

    // Wait until the next event or the remaining timeout
    {
      MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);

      if (should_terminate() || _collector->_full_gc_requested) {
        return;
      }
      set_gogc_flag(gogc_gogc_wants_token);   // to provoke notifies
      assert(t_millis == 0 || wait_time_millis > 0, "Sanity");
      CGC_lock->wait_without_safepoint_check(wait_time_millis);
      clear_gogc_flag(gogc_gogc_wants_token);
      assert(!gogc_flag_is_set(gogc_gogc_has_token | gogc_gogc_wants_token),
             "Should not be set");
    }

    if(t_millis != 0 && os::elapsedTime() >= end_time_secs) {
      break;
    }

    unsigned int after_count;
    {
      MutexLocker hl(Heap_lock, Mutex::_no_safepoint_check_flag);
      after_count = heap->total_collections();
    }

    if(before_count != after_count) {
      break;
    }

    if(++loop_count == 0) {
      log_warning(gc)("wait_on_gogc_lock_for_scavenge() has looped %u times", loop_count - 1);
    }
  }
}

void goConcurrentMarkSweepThread::sleepBeforeNextCycle() {
  while (!should_terminate()) {
    if(gogcWaitDuration >= 0) {
      wait_on_gogc_lock_for_scavenge(gogcWaitDuration);
    } else {
      wait_on_gogc_lock(gogcCheckInterval);
    }
    if (_collector->shouldConcurrentCollect()) {
      return;
    }
  }
}
