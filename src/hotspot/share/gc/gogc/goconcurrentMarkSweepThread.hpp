#ifndef SHARE_GC_gogc_goConcurrentMarkSweepThread_HPP
#define SHARE_GC_gogc_goConcurrentMarkSweepThread_HPP

#include "gc/gogc/concurrentMarkSweepGeneration.hpp"
#include "gc/shared/concurrentGCThread.hpp"
#include "runtime/thread.hpp"

class ConcurrentMarkSweepGeneration;
class gogcCollector;

// The Concurrent Mark Sweep GC Thread
class goConcurrentMarkSweepThread: public ConcurrentGCThread {
  friend class VMStructs;
  friend class ConcurrentMarkSweepGeneration;  
  friend class gogcCollector;

 private:
  static goConcurrentMarkSweepThread* _gogct;
  static gogcCollector*              _collector;

  enum gogc_flag_type {
    gogc_nil             = NoBits,
    gogc_gogc_wants_token = nth_bit(0),
    gogc_gogc_has_token   = nth_bit(1),
    gogc_vm_wants_token  = nth_bit(2),
    gogc_vm_has_token    = nth_bit(3)
  };

  static int _gogc_flag;

  static bool gogc_flag_is_set(int b)        { return (_gogc_flag & b) != 0;   }
  static bool set_gogc_flag(int b)           { return (_gogc_flag |= b) != 0;  }
  static bool clear_gogc_flag(int b)         { return (_gogc_flag &= ~b) != 0; }
  void sleepBeforeNextCycle();

  // gogc thread should yield for a young gen collection and direct allocations
  static char _pad_1[64 - sizeof(jint)];    // prevent cache-line sharing
  static volatile jint _pending_yields;
  static char _pad_2[64 - sizeof(jint)];    // prevent cache-line sharing

  // debugging
  void verify_ok_to_terminate() const PRODUCT_RETURN;

  void run_service();
  void stop_service();

 public:
  // Constructor
  goConcurrentMarkSweepThread(gogcCollector* collector);

  static void threads_do(ThreadClosure* tc);

  // Printing
  static void print_all_on(outputStream* st);
  static void print_all()                             { print_all_on(tty); }

  // Returns the gogc Thread
  static goConcurrentMarkSweepThread* gogct()    { return _gogct; }
  static gogcCollector*         collector()    { return _collector;  }

  // Create and start the gogc Thread, or stop it on shutdown
  static goConcurrentMarkSweepThread* start(gogcCollector* collector);

  // Synchronization using gogc token
  static void synchronize(bool is_gogc_thread);
  static void desynchronize(bool is_gogc_thread);
  static bool vm_thread_has_gogc_token() {
    return gogc_flag_is_set(gogc_vm_has_token);
  }
  static bool gogc_thread_has_gogc_token() {
    return gogc_flag_is_set(gogc_gogc_has_token);
  }
  static bool vm_thread_wants_gogc_token() {
    return gogc_flag_is_set(gogc_vm_wants_token);
  }
  static bool gogc_thread_wants_gogc_token() {
    return gogc_flag_is_set(gogc_gogc_wants_token);
  }

  // Wait on gogc lock until the next synchronous GC
  // or given timeout, whichever is earlier. A timeout value
  // of 0 indicates that there is no upper bound on the wait time.
  // A concurrent full gc request terminates the wait.
  void wait_on_gogc_lock(long t_millis);

  // Wait on gogc lock until the next synchronous GC
  // or given timeout, whichever is earlier. A timeout value
  // of 0 indicates that there is no upper bound on the wait time.
  // A concurrent full gc request terminates the wait.
  void wait_on_gogc_lock_for_scavenge(long t_millis);

  // The gogc thread will yield during the work portion of its cycle
  // only when requested to.
  // A synchronous request is used for young gen collections and
  // for direct allocations.  The requesting thread increments
  // _pending_yields at the beginning of an operation, and decrements
  // _pending_yields when that operation is completed.
  // In turn, the gogc thread yields when _pending_yields is positive,
  // and continues to yield until the value reverts to 0.

  static void increment_pending_yields()   {
    Atomic::inc(&_pending_yields);
    assert(_pending_yields >= 0, "can't be negative");
  }
  static void decrement_pending_yields()   {
    Atomic::dec(&_pending_yields);
    assert(_pending_yields >= 0, "can't be negative");
  }
  static bool should_yield()   { return _pending_yields > 0; }
};

// For scoped increment/decrement of (synchronous) yield requests
class gogcSynchronousYieldRequest: public StackObj {
 public:
  gogcSynchronousYieldRequest() {
    goConcurrentMarkSweepThread::increment_pending_yields();
  }
  ~gogcSynchronousYieldRequest() {
    goConcurrentMarkSweepThread::decrement_pending_yields();
  }
};

// Used to emit a warning in case of unexpectedly excessive
// looping (in "apparently endless loops") in gogc code.
class gogcLoopCountWarn: public StackObj {
 private:
  const char* _src;
  const char* _msg;
  const intx  _threshold;
  intx        _ticks;

 public:
  inline gogcLoopCountWarn(const char* src, const char* msg,
                          const intx threshold) :
    _src(src), _msg(msg), _threshold(threshold), _ticks(0) { }

  inline void tick() {
    _ticks++;
    if (gogcLoopWarn && _ticks % _threshold == 0) {
      log_warning(gc)("%s has looped " INTX_FORMAT " times %s", _src, _ticks, _msg);
    }
  }
};

#endif // SHARE_GC_gogc_goConcurrentMarkSweepThread_HPP
