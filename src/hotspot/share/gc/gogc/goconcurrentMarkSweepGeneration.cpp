#include "precompiled.hpp"
#include "classfile/classLoaderDataGraph.hpp"
#include "classfile/systemDictionary.hpp"
#include "code/codeCache.hpp"
#include "gc/gogc/gogcGCStats.hpp"
#include "gc/gogc/gogcHeap.hpp"
#include "gc/gogc/gogcOopClosures.inline.hpp"
#include "gc/gogc/gogcVMOperations.hpp"
#include "gc/gogc/compactibleFreeListSpace.hpp"
#include "gc/gogc/goConcurrentMarkSweepGeneration.inline.hpp"
#include "gc/gogc/concurrentMarkSweepThread.hpp"
#include "gc/gogc/parNewGeneration.hpp"
#include "gc/gogc/promotionInfo.inline.hpp"
#include "gc/serial/genMarkSweep.hpp"
#include "gc/serial/tenuredGeneration.hpp"
#include "gc/shared/adaptiveSizePolicy.hpp"
#include "gc/shared/cardGeneration.inline.hpp"
#include "gc/shared/cardTableRS.hpp"
#include "gc/shared/collectedHeap.inline.hpp"
#include "gc/shared/collectorCounters.hpp"
#include "gc/shared/gcLocker.hpp"
#include "gc/shared/gcPolicyCounters.hpp"
#include "gc/shared/gcTimer.hpp"
#include "gc/shared/gcTrace.hpp"
#include "gc/shared/gcTraceTime.inline.hpp"
#include "gc/shared/genCollectedHeap.hpp"
#include "gc/shared/genOopClosures.inline.hpp"
#include "gc/shared/isGCActiveMark.hpp"
#include "gc/shared/owstTaskTerminator.hpp"
#include "gc/shared/referencePolicy.hpp"
#include "gc/shared/referenceProcessorPhaseTimes.hpp"
#include "gc/shared/space.inline.hpp"
#include "gc/shared/strongRootsScope.hpp"
#include "gc/shared/taskqueue.inline.hpp"
#include "gc/shared/weakProcessor.hpp"
#include "gc/shared/workerPolicy.hpp"
#include "logging/log.hpp"
#include "logging/logStream.hpp"
#include "memory/allocation.hpp"
#include "memory/binaryTreeDictionary.inline.hpp"
#include "memory/iterator.inline.hpp"
#include "memory/padded.hpp"
#include "memory/resourceArea.hpp"
#include "memory/universe.hpp"
#include "oops/access.inline.hpp"
#include "oops/oop.inline.hpp"
#include "prims/jvmtiExport.hpp"
#include "runtime/atomic.hpp"
#include "runtime/flags/flagSetting.hpp"
#include "runtime/globals_extension.hpp"
#include "runtime/handles.inline.hpp"
#include "runtime/java.hpp"
#include "runtime/orderAccess.hpp"
#include "runtime/timer.hpp"
#include "runtime/vmThread.hpp"
#include "services/memoryService.hpp"
#include "services/runtimeService.hpp"
#include "utilities/align.hpp"
#include "utilities/stack.inline.hpp"
#if INCLUDE_JVMCI
#include "jvmci/jvmci.hpp"
#endif

// statics
gogcCollector* goConcurrentMarkSweepGeneration::_collector = NULL;
bool gogcCollector::_full_gc_requested = false;
GCCause::Cause gogcCollector::_full_gc_cause = GCCause::_no_gc;

class gogcTokenSync: public StackObj {
 private:
  bool _is_gogc_thread;
 public:
  gogcTokenSync(bool is_gogc_thread):
    _is_gogc_thread(is_gogc_thread) {
    assert(is_gogc_thread == Thread::current()->is_ConcurrentGC_thread(),
           "Incorrect argument to constructor");
    ConcurrentMarkSweepThread::synchronize(_is_gogc_thread);
  }

  ~gogcTokenSync() {
    assert(_is_gogc_thread ?
             ConcurrentMarkSweepThread::gogc_thread_has_gogc_token() :
             ConcurrentMarkSweepThread::vm_thread_has_gogc_token(),
          "Incorrect state");
    ConcurrentMarkSweepThread::desynchronize(_is_gogc_thread);
  }
};

// Convenience class that does a gogcTokenSync, and then acquires
// upto three locks.
class gogcTokenSyncWithLocks: public gogcTokenSync {
 private:
  // Note: locks are acquired in textual declaration order
  // and released in the opposite order
  MutexLocker _locker1, _locker2, _locker3;
 public:
  gogcTokenSyncWithLocks(bool is_gogc_thread, Mutex* mutex1,
                        Mutex* mutex2 = NULL, Mutex* mutex3 = NULL):
    gogcTokenSync(is_gogc_thread),
    _locker1(mutex1, Mutex::_no_safepoint_check_flag),
    _locker2(mutex2, Mutex::_no_safepoint_check_flag),
    _locker3(mutex3, Mutex::_no_safepoint_check_flag)
  { }
};


NOT_PRODUCT(CompactibleFreeListSpace* debug_gogc_space;)

// This struct contains per-thread things necessary to support parallel
// young-gen collection.
class gogcParGCThreadState: public CHeapObj<mtGC> {
 public:
  CompactibleFreeListSpaceLAB lab;
  PromotionInfo promo;

  // Constructor.
  gogcParGCThreadState(CompactibleFreeListSpace* cfls) : lab(cfls) {
    promo.setSpace(cfls);
  }
};

goConcurrentMarkSweepGeneration::goConcurrentMarkSweepGeneration(
     ReservedSpace rs,
     size_t initial_byte_size,
     size_t min_byte_size,
     size_t max_byte_size,
     CardTableRS* ct) :
  CardGeneration(rs, initial_byte_size, ct),
  _dilatation_factor(((double)MinChunkSize)/((double)(CollectedHeap::min_fill_size()))),
  _did_compact(false)
{
  HeapWord* bottom = (HeapWord*) _virtual_space.low();
  HeapWord* end    = (HeapWord*) _virtual_space.high();

  _direct_allocated_words = 0;
  NOT_PRODUCT(
    _numObjectsPromoted = 0;
    _numWordsPromoted = 0;
    _numObjectsAllocated = 0;
    _numWordsAllocated = 0;
  )

  _gogcSpace = new CompactibleFreeListSpace(_bts, MemRegion(bottom, end));
  NOT_PRODUCT(debug_gogc_space = _gogcSpace;)
  _gogcSpace->_old_gen = this;

  _gc_stats = new gogcGCStats();

  // Verify the assumption that FreeChunk::_prev and OopDesc::_klass
  // offsets match. The ability to tell free chunks from objects
  // depends on this property.
  debug_only(
    FreeChunk* junk = NULL;
    assert(UseCompressedClassPointers ||
           junk->prev_addr() == (void*)(oop(junk)->klass_addr()),
           "Offset of FreeChunk::_prev within FreeChunk must match"
           "  that of OopDesc::_klass within OopDesc");
  )

  _par_gc_thread_states = NEW_C_HEAP_ARRAY(gogcParGCThreadState*, ParallelGCThreads, mtGC);
  for (uint i = 0; i < ParallelGCThreads; i++) {
    _par_gc_thread_states[i] = new gogcParGCThreadState(gogcSpace());
  }

  _incremental_collection_failed = false;

  assert(MinChunkSize >= CollectedHeap::min_fill_size(), "just checking");
  assert(_dilatation_factor >= 1.0, "from previous assert");

  initialize_performance_counters(min_byte_size, max_byte_size);
}


void goConcurrentMarkSweepGeneration::init_initiating_occupancy(intx io, uintx tr) {
  assert(io <= 100 && tr <= 100, "Check the arguments");
  if (io >= 0) {
    _initiating_occupancy = (double)io / 100.0;
  } else {
    _initiating_occupancy = ((100 - MinHeapFreeRatio) +
                             (double)(tr * MinHeapFreeRatio) / 100.0)
                            / 100.0;
  }
}

void goConcurrentMarkSweepGeneration::ref_processor_init() {
  assert(collector() != NULL, "no collector");
  collector()->ref_processor_init();
}

void gogcCollector::ref_processor_init() {
  if (_ref_processor == NULL) {
    // Allocate and initialize a reference processor
    _ref_processor =
      new ReferenceProcessor(&_span_based_discoverer,
                             (ParallelGCThreads > 1) && ParallelRefProcEnabled, // mt processing
                             ParallelGCThreads,                      // mt processing degree
                             _gogcGen->refs_discovery_is_mt(),        // mt discovery
                             MAX2(ConcGCThreads, ParallelGCThreads), // mt discovery degree
                             _gogcGen->refs_discovery_is_atomic(),    // discovery is not atomic
                             &_is_alive_closure,                     // closure for liveness info
                             false);                                 // disable adjusting number of processing threads
    // Initialize the _ref_processor field of gogcGen
    _gogcGen->set_ref_processor(_ref_processor);

  }
}

AdaptiveSizePolicy* gogcCollector::size_policy() {
  return gogcHeap::heap()->size_policy();
}

void goConcurrentMarkSweepGeneration::initialize_performance_counters(size_t min_old_size,
                                                                    size_t max_old_size) {

  const char* gen_name = "old";
  // Generation Counters - generation 1, 1 subspace
  _gen_counters = new GenerationCounters(gen_name, 1, 1,
      min_old_size, max_old_size, &_virtual_space);

  _space_counters = new GSpaceCounters(gen_name, 0,
                                       _virtual_space.reserved_size(),
                                       this, _gen_counters);
}

gogcStats::gogcStats(goConcurrentMarkSweepGeneration* gogc_gen, unsigned int alpha):
  _gogc_gen(gogc_gen)
{
  assert(alpha <= 100, "bad value");
  _saved_alpha = alpha;

  // Initialize the alphas to the bootstrap value of 100.
  _gc0_alpha = _gogc_alpha = 100;

  _gogc_begin_time.update();
  _gogc_end_time.update();

  _gc0_duration = 0.0;
  _gc0_period = 0.0;
  _gc0_promoted = 0;

  _gogc_duration = 0.0;
  _gogc_period = 0.0;
  _gogc_allocated = 0;

  _gogc_used_at_gc0_begin = 0;
  _gogc_used_at_gc0_end = 0;
  _allow_duty_cycle_reduction = false;
  _valid_bits = 0;
}

double gogcStats::gogc_free_adjustment_factor(size_t free) const {
  // TBD: CR 6909490
  return 1.0;
}

void gogcStats::adjust_gogc_free_adjustment_factor(bool fail, size_t free) {
}

// If promotion failure handling is on use
// the padded average size of the promotion for each
// young generation collection.
double gogcStats::time_until_gogc_gen_full() const {
  size_t gogc_free = _gogc_gen->gogcSpace()->free();
  gogcHeap* heap = gogcHeap::heap();
  size_t expected_promotion = MIN2(heap->young_gen()->capacity(),
                                   (size_t) _gogc_gen->gc_stats()->avg_promoted()->padded_average());
  if (gogc_free > expected_promotion) {
    // Start a gogc collection if there isn't enough space to promote
    // for the next young collection.  Use the padded average as
    // a safety factor.
    gogc_free -= expected_promotion;

    // Adjust by the safety factor.
    double gogc_free_dbl = (double)gogc_free;
    double gogc_adjustment = (100.0 - gogcIncrementalSafetyFactor) / 100.0;
    // Apply a further correction factor which tries to adjust
    // for recent occurance of concurrent mode failures.
    gogc_adjustment = gogc_adjustment * gogc_free_adjustment_factor(gogc_free);
    gogc_free_dbl = gogc_free_dbl * gogc_adjustment;

    log_trace(gc)("gogcStats::time_until_gogc_gen_full: gogc_free " SIZE_FORMAT " expected_promotion " SIZE_FORMAT,
                  gogc_free, expected_promotion);
    log_trace(gc)("  gogc_free_dbl %f gogc_consumption_rate %f", gogc_free_dbl, gogc_consumption_rate() + 1.0);
    // Add 1 in case the consumption rate goes to zero.
    return gogc_free_dbl / (gogc_consumption_rate() + 1.0);
  }
  return 0.0;
}

// Compare the duration of the gogc collection to the
// time remaining before the gogc generation is empty.
// Note that the time from the start of the gogc collection
// to the start of the gogc sweep (less than the total
// duration of the gogc collection) can be used.  This
// has been tried and some applications experienced
// promotion failures early in execution.  This was
// possibly because the averages were not accurate
// enough at the beginning.
double gogcStats::time_until_gogc_start() const {
  // We add "gc0_period" to the "work" calculation
  // below because this query is done (mostly) at the
  // end of a scavenge, so we need to conservatively
  // account for that much possible delay
  // in the query so as to avoid concurrent mode failures
  // due to starting the collection just a wee bit too
  // late.
  double work = gogc_duration() + gc0_period();
  double deadline = time_until_gogc_gen_full();
  // If a concurrent mode failure occurred recently, we want to be
  // more conservative and halve our expected time_until_gogc_gen_full()
  if (work > deadline) {
    log_develop_trace(gc)("gogcCollector: collect because of anticipated promotion before full %3.7f + %3.7f > %3.7f ",
                          gogc_duration(), gc0_period(), time_until_gogc_gen_full());
    return 0.0;
  }
  return work - deadline;
}

#ifndef PRODUCT
void gogcStats::print_on(outputStream *st) const {
  st->print(" gc0_alpha=%d,gogc_alpha=%d", _gc0_alpha, _gogc_alpha);
  st->print(",gc0_dur=%g,gc0_per=%g,gc0_promo=" SIZE_FORMAT,
               gc0_duration(), gc0_period(), gc0_promoted());
  st->print(",gogc_dur=%g,gogc_per=%g,gogc_alloc=" SIZE_FORMAT,
            gogc_duration(), gogc_period(), gogc_allocated());
  st->print(",gogc_since_beg=%g,gogc_since_end=%g",
            gogc_time_since_begin(), gogc_time_since_end());
  st->print(",gogc_used_beg=" SIZE_FORMAT ",gogc_used_end=" SIZE_FORMAT,
            _gogc_used_at_gc0_begin, _gogc_used_at_gc0_end);

  if (valid()) {
    st->print(",promo_rate=%g,gogc_alloc_rate=%g",
              promotion_rate(), gogc_allocation_rate());
    st->print(",gogc_consumption_rate=%g,time_until_full=%g",
              gogc_consumption_rate(), time_until_gogc_gen_full());
  }
  st->cr();
}
#endif // #ifndef PRODUCT

gogcCollector::CollectorState gogcCollector::_collectorState =
                             gogcCollector::Idling;
bool gogcCollector::_foregroundGCIsActive = false;
bool gogcCollector::_foregroundGCShouldWait = false;

gogcCollector::gogcCollector(goConcurrentMarkSweepGeneration* gogcGen,
                           CardTableRS*                   ct):
  _overflow_list(NULL),
  _conc_workers(NULL),     // may be set later
  _completed_initialization(false),
  _collection_count_start(0),
  _should_unload_classes(gogcClassUnloadingEnabled),
  _concurrent_cycles_since_last_unload(0),
  _roots_scanning_options(GenCollectedHeap::SO_None),
  _verification_mark_bm(0, Mutex::leaf + 1, "gogc_verification_mark_bm_lock"),
  _verifying(false),
  _inter_sweep_estimate(gogc_SweepWeight, gogc_SweepPadding),
  _intra_sweep_estimate(gogc_SweepWeight, gogc_SweepPadding),
  _gc_tracer_cm(new (ResourceObj::C_HEAP, mtGC) gogcTracer()),
  _gc_timer_cm(new (ResourceObj::C_HEAP, mtGC) ConcurrentGCTimer()),
  _gogc_start_registered(false),
  _gogcGen(gogcGen),
  // Adjust span to cover old (gogc) gen
  _span(gogcGen->reserved()),
  _ct(ct),
  _markBitMap(0, Mutex::leaf + 1, "gogc_markBitMap_lock"),
  _modUnionTable((CardTable::card_shift - LogHeapWordSize),
                 -1 /* lock-free */, "No_lock" /* dummy */),
  _restart_addr(NULL),
  _ser_pmc_preclean_ovflw(0),
  _ser_pmc_remark_ovflw(0),
  _par_pmc_remark_ovflw(0),
  _ser_kac_preclean_ovflw(0),
  _ser_kac_ovflw(0),
  _par_kac_ovflw(0),
#ifndef PRODUCT
  _num_par_pushes(0),
#endif
  _span_based_discoverer(_span),
  _ref_processor(NULL),    // will be set later
  // Construct the is_alive_closure with _span & markBitMap
  _is_alive_closure(_span, &_markBitMap),
  _modUnionClosurePar(&_modUnionTable),
  _between_prologue_and_epilogue(false),
  _abort_preclean(false),
  _start_sampling(false),
  _stats(gogcGen),
  _eden_chunk_lock(new Mutex(Mutex::leaf + 1, "gogc_eden_chunk_lock", true,
                             //verify that this lock should be acquired with safepoint check.
                             Monitor::_safepoint_check_never)),
  _eden_chunk_array(NULL),     // may be set in ctor body
  _eden_chunk_index(0),        // -- ditto --
  _eden_chunk_capacity(0),     // -- ditto --
  _survivor_chunk_array(NULL), // -- ditto --
  _survivor_chunk_index(0),    // -- ditto --
  _survivor_chunk_capacity(0), // -- ditto --
  _survivor_plab_array(NULL)   // -- ditto --
{
  // Now expand the span and allocate the collection support structures
  // (MUT, marking bit map etc.) to cover both generations subject to
  // collection.

  // For use by dirty card to oop closures.
  _gogcGen->gogcSpace()->set_collector(this);

  // Allocate MUT and marking bit map
  {
    MutexLocker x(_markBitMap.lock(), Mutex::_no_safepoint_check_flag);
    if (!_markBitMap.allocate(_span)) {
      log_warning(gc)("Failed to allocate gogc Bit Map");
      return;
    }
    assert(_markBitMap.covers(_span), "_markBitMap inconsistency?");
  }
  {
    _modUnionTable.allocate(_span);
    assert(_modUnionTable.covers(_span), "_modUnionTable inconsistency?");
  }

  if (!_markStack.allocate(MarkStackSize)) {
    log_warning(gc)("Failed to allocate gogc Marking Stack");
    return;
  }

  // Support for multi-threaded concurrent phases
  if (gogcConcurrentMTEnabled) {
    if (FLAG_IS_DEFAULT(ConcGCThreads)) {
      // just for now
      FLAG_SET_DEFAULT(ConcGCThreads, (ParallelGCThreads + 3) / 4);
    }
    if (ConcGCThreads > 1) {
      _conc_workers = new YieldingFlexibleWorkGang("gogc Thread",
                                 ConcGCThreads, true);
      if (_conc_workers == NULL) {
        log_warning(gc)("GC/gogc: _conc_workers allocation failure: forcing -gogcConcurrentMTEnabled");
        gogcConcurrentMTEnabled = false;
      } else {
        _conc_workers->initialize_workers();
      }
    } else {
      gogcConcurrentMTEnabled = false;
    }
  }
  if (!gogcConcurrentMTEnabled) {
    ConcGCThreads = 0;
  } else {
    // Turn off gogcCleanOnEnter optimization temporarily for
    // the MT case where it's not fixed yet; see 6178663.
    gogcCleanOnEnter = false;
  }
  assert((_conc_workers != NULL) == (ConcGCThreads > 1),
         "Inconsistency");
  log_debug(gc)("ConcGCThreads: %u", ConcGCThreads);
  log_debug(gc)("ParallelGCThreads: %u", ParallelGCThreads);

  // Parallel task queues; these are shared for the
  // concurrent and stop-world phases of gogc, but
  // are not shared with parallel scavenge (ParNew).
  {
    uint i;
    uint num_queues = MAX2(ParallelGCThreads, ConcGCThreads);

    if ((gogcParallelRemarkEnabled || gogcConcurrentMTEnabled
         || ParallelRefProcEnabled)
        && num_queues > 0) {
      _task_queues = new OopTaskQueueSet(num_queues);
      if (_task_queues == NULL) {
        log_warning(gc)("task_queues allocation failure.");
        return;
      }
      typedef Padded<OopTaskQueue> PaddedOopTaskQueue;
      for (i = 0; i < num_queues; i++) {
        PaddedOopTaskQueue *q = new PaddedOopTaskQueue();
        if (q == NULL) {
          log_warning(gc)("work_queue allocation failure.");
          return;
        }
        _task_queues->register_queue(i, q);
      }
      for (i = 0; i < num_queues; i++) {
        _task_queues->queue(i)->initialize();
      }
    }
  }

  _gogcGen ->init_initiating_occupancy(gogcInitiatingOccupancyFraction, gogcTriggerRatio);

  // Clip gogcBootstrapOccupancy between 0 and 100.
  _bootstrap_occupancy = gogcBootstrapOccupancy / 100.0;

  // Now tell gogc generations the identity of their collector
  goConcurrentMarkSweepGeneration::set_collector(this);

  // Create & start a gogc thread for this gogc collector
  _gogcThread = ConcurrentMarkSweepThread::start(this);
  assert(gogcThread() != NULL, "gogc Thread should have been created");
  assert(gogcThread()->collector() == this,
         "gogc Thread should refer to this gen");
  assert(CGC_lock != NULL, "Where's the CGC_lock?");

  // Support for parallelizing young gen rescan
  gogcHeap* heap = gogcHeap::heap();
  _young_gen = heap->young_gen();
  if (heap->supports_inline_contig_alloc()) {
    _top_addr = heap->top_addr();
    _end_addr = heap->end_addr();
    assert(_young_gen != NULL, "no _young_gen");
    _eden_chunk_index = 0;
    _eden_chunk_capacity = (_young_gen->max_capacity() + gogcSamplingGrain) / gogcSamplingGrain;
    _eden_chunk_array = NEW_C_HEAP_ARRAY(HeapWord*, _eden_chunk_capacity, mtGC);
  }

  // Support for parallelizing survivor space rescan
  if ((gogcParallelRemarkEnabled && gogcParallelSurvivorRemarkEnabled) || gogcParallelInitialMarkEnabled) {
    const size_t max_plab_samples =
      _young_gen->max_survivor_size() / (PLAB::min_size() * HeapWordSize);

    _survivor_plab_array  = NEW_C_HEAP_ARRAY(ChunkArray, ParallelGCThreads, mtGC);
    _survivor_chunk_array = NEW_C_HEAP_ARRAY(HeapWord*, max_plab_samples, mtGC);
    _cursor               = NEW_C_HEAP_ARRAY(size_t, ParallelGCThreads, mtGC);
    _survivor_chunk_capacity = max_plab_samples;
    for (uint i = 0; i < ParallelGCThreads; i++) {
      HeapWord** vec = NEW_C_HEAP_ARRAY(HeapWord*, max_plab_samples, mtGC);
      ChunkArray* cur = ::new (&_survivor_plab_array[i]) ChunkArray(vec, max_plab_samples);
      assert(cur->end() == 0, "Should be 0");
      assert(cur->array() == vec, "Should be vec");
      assert(cur->capacity() == max_plab_samples, "Error");
    }
  }

  NOT_PRODUCT(_overflow_counter = gogcMarkStackOverflowInterval;)
  _gc_counters = new CollectorCounters("gogc full collection pauses", 1);
  _cgc_counters = new CollectorCounters("gogc concurrent cycle pauses", 2);
  _completed_initialization = true;
  _inter_sweep_timer.start();  // start of time
}

const char* goConcurrentMarkSweepGeneration::name() const {
  return "concurrent mark-sweep generation";
}
void goConcurrentMarkSweepGeneration::update_counters() {
  if (UsePerfData) {
    _space_counters->update_all();
    _gen_counters->update_all();
  }
}

// this is an optimized version of update_counters(). it takes the
// used value as a parameter rather than computing it.
//
void goConcurrentMarkSweepGeneration::update_counters(size_t used) {
  if (UsePerfData) {
    _space_counters->update_used(used);
    _space_counters->update_capacity();
    _gen_counters->update_all();
  }
}

void goConcurrentMarkSweepGeneration::print() const {
  Generation::print();
  gogcSpace()->print();
}

#ifndef PRODUCT
void goConcurrentMarkSweepGeneration::print_statistics() {
  gogcSpace()->printFLCensus(0);
}
#endif

size_t
goConcurrentMarkSweepGeneration::contiguous_available() const {

  return MAX2(_virtual_space.uncommitted_size(), unsafe_max_alloc_nogc());
}

size_t
goConcurrentMarkSweepGeneration::unsafe_max_alloc_nogc() const {
  return _gogcSpace->max_alloc_in_words() * HeapWordSize;
}

size_t goConcurrentMarkSweepGeneration::max_available() const {
  return free() + _virtual_space.uncommitted_size();
}

bool goConcurrentMarkSweepGeneration::promotion_attempt_is_safe(size_t max_promotion_in_bytes) const {
  size_t available = max_available();
  size_t av_promo  = (size_t)gc_stats()->avg_promoted()->padded_average();
  bool   res = (available >= av_promo) || (available >= max_promotion_in_bytes);
  log_trace(gc, promotion)("gogc: promo attempt is%s safe: available(" SIZE_FORMAT ") %s av_promo(" SIZE_FORMAT "), max_promo(" SIZE_FORMAT ")",
                           res? "":" not", available, res? ">=":"<", av_promo, max_promotion_in_bytes);
  return res;
}

// At a promotion failure dump information on block layout in heap
// (gogc old generation).
void goConcurrentMarkSweepGeneration::promotion_failure_occurred() {
  Log(gc, promotion) log;
  if (log.is_trace()) {
    LogStream ls(log.trace());
    gogcSpace()->dump_at_safepoint_with_locks(collector(), &ls);
  }
}

void goConcurrentMarkSweepGeneration::reset_after_compaction() {

  for (uint i = 0; i < ParallelGCThreads; i++) {
    _par_gc_thread_states[i]->promo.reset();
  }
}

void goConcurrentMarkSweepGeneration::compute_new_size() {
  assert_locked_or_safepoint(Heap_lock);

  if (incremental_collection_failed()) {
    clear_incremental_collection_failed();
    grow_to_reserved();
    return;
  }


  CardGeneration::compute_new_size();

  if (did_compact()) {
    gogcSpace()->reset_after_compaction();
  }
}

void goConcurrentMarkSweepGeneration::compute_new_size_free_list() {
  assert_locked_or_safepoint(Heap_lock);

  if (incremental_collection_failed()) {
    clear_incremental_collection_failed();
    grow_to_reserved();
    return;
  }

  double free_percentage = ((double) free()) / capacity();
  double desired_free_percentage = (double) MinHeapFreeRatio / 100;
  double maximum_free_percentage = (double) MaxHeapFreeRatio / 100;

  if (free_percentage < desired_free_percentage) {
    size_t desired_capacity = (size_t)(used() / ((double) 1 - desired_free_percentage));
    assert(desired_capacity >= capacity(), "invalid expansion size");
    size_t expand_bytes = MAX2(desired_capacity - capacity(), MinHeapDeltaBytes);
    Log(gc) log;
    if (log.is_trace()) {
      size_t desired_capacity = (size_t)(used() / ((double) 1 - desired_free_percentage));
      log.trace("From compute_new_size: ");
      log.trace("  Free fraction %f", free_percentage);
      log.trace("  Desired free fraction %f", desired_free_percentage);
      log.trace("  Maximum free fraction %f", maximum_free_percentage);
      log.trace("  Capacity " SIZE_FORMAT, capacity() / 1000);
      log.trace("  Desired capacity " SIZE_FORMAT, desired_capacity / 1000);
      gogcHeap* heap = gogcHeap::heap();
      size_t young_size = heap->young_gen()->capacity();
      log.trace("  Young gen size " SIZE_FORMAT, young_size / 1000);
      log.trace("  unsafe_max_alloc_nogc " SIZE_FORMAT, unsafe_max_alloc_nogc() / 1000);
      log.trace("  contiguous available " SIZE_FORMAT, contiguous_available() / 1000);
      log.trace("  Expand by " SIZE_FORMAT " (bytes)", expand_bytes);
    }
    // safe if expansion fails
    expand_for_gc_cause(expand_bytes, 0, gogcExpansionCause::_satisfy_free_ratio);
    log.trace("  Expanded free fraction %f", ((double) free()) / capacity());
  } else {
    size_t desired_capacity = (size_t)(used() / ((double) 1 - desired_free_percentage));
    assert(desired_capacity <= capacity(), "invalid expansion size");
    size_t shrink_bytes = capacity() - desired_capacity;
    // Don't shrink unless the delta is greater than the minimum shrink we want
    if (shrink_bytes >= MinHeapDeltaBytes) {
      shrink_free_list_by(shrink_bytes);
    }
  }
}

Mutex* goConcurrentMarkSweepGeneration::freelistLock() const {
  return gogcSpace()->freelistLock();
}

HeapWord* goConcurrentMarkSweepGeneration::allocate(size_t size, bool tlab) {
  gogcSynchronousYieldRequest yr;
  MutexLocker x(freelistLock(), Mutex::_no_safepoint_check_flag);
  return have_lock_and_allocate(size, tlab);
}

HeapWord* goConcurrentMarkSweepGeneration::have_lock_and_allocate(size_t size,
                                                                bool   tlab /* ignored */) {
  assert_lock_strong(freelistLock());
  size_t adjustedSize = CompactibleFreeListSpace::adjustObjectSize(size);
  HeapWord* res = gogcSpace()->allocate(adjustedSize);
  // Allocate the object live (grey) if the background collector has
  // started marking. This is necessary because the marker may
  // have passed this address and consequently this object will
  // not otherwise be greyed and would be incorrectly swept up.
  // Note that if this object contains references, the writing
  // of those references will dirty the card containing this object
  // allowing the object to be blackened (and its references scanned)
  // either during a preclean phase or at the final checkpoint.
  if (res != NULL) {
    // We may block here with an uninitialized object with
    // its mark-bit or P-bits not yet set. Such objects need
    // to be safely navigable by block_start().
    assert(oop(res)->klass_or_null() == NULL, "Object should be uninitialized here.");
    assert(!((FreeChunk*)res)->is_free(), "Error, block will look free but show wrong size");
    collector()->direct_allocated(res, adjustedSize);
    _direct_allocated_words += adjustedSize;
    // allocation counters
    NOT_PRODUCT(
      _numObjectsAllocated++;
      _numWordsAllocated += (int)adjustedSize;
    )
  }
  return res;
}

// In the case of direct allocation by mutators in a generation that
// is being concurrently collected, the object must be allocated
// live (grey) if the background collector has started marking.
// This is necessary because the marker may
// have passed this address and consequently this object will
// not otherwise be greyed and would be incorrectly swept up.
// Note that if this object contains references, the writing
// of those references will dirty the card containing this object
// allowing the object to be blackened (and its references scanned)
// either during a preclean phase or at the final checkpoint.
void gogcCollector::direct_allocated(HeapWord* start, size_t size) {
  assert(_markBitMap.covers(start, size), "Out of bounds");
  if (_collectorState >= Marking) {
    MutexLocker y(_markBitMap.lock(),
                  Mutex::_no_safepoint_check_flag);
    // [see comments preceding SweepClosure::do_blk() below for details]
    //
    // Can the P-bits be deleted now?  JJJ
    //
    // 1. need to mark the object as live so it isn't collected
    // 2. need to mark the 2nd bit to indicate the object may be uninitialized
    // 3. need to mark the end of the object so marking, precleaning or sweeping
    //    can skip over uninitialized or unparsable objects. An allocated
    //    object is considered uninitialized for our purposes as long as
    //    its klass word is NULL.  All old gen objects are parsable
    //    as soon as they are initialized.)
    _markBitMap.mark(start);          // object is live
    _markBitMap.mark(start + 1);      // object is potentially uninitialized?
    _markBitMap.mark(start + size - 1);
                                      // mark end of object
  }
  // check that oop looks uninitialized
  assert(oop(start)->klass_or_null() == NULL, "_klass should be NULL");
}

void gogcCollector::promoted(bool par, HeapWord* start,
                            bool is_obj_array, size_t obj_size) {
  assert(_markBitMap.covers(start), "Out of bounds");
  // See comment in direct_allocated() about when objects should
  // be allocated live.
  if (_collectorState >= Marking) {
    // we already hold the marking bit map lock, taken in
    // the prologue
    if (par) {
      _markBitMap.par_mark(start);
    } else {
      _markBitMap.mark(start);
    }
    // We don't need to mark the object as uninitialized (as
    // in direct_allocated above) because this is being done with the
    // world stopped and the object will be initialized by the
    // time the marking, precleaning or sweeping get to look at it.
    // But see the code for copying objects into the gogc generation,
    // where we need to ensure that concurrent readers of the
    // block offset table are able to safely navigate a block that
    // is in flux from being free to being allocated (and in
    // transition while being copied into) and subsequently
    // becoming a bona-fide object when the copy/promotion is complete.
    assert(SafepointSynchronize::is_at_safepoint(),
           "expect promotion only at safepoints");

    if (_collectorState < Sweeping) {
      // Mark the appropriate cards in the modUnionTable, so that
      // this object gets scanned before the sweep. If this is
      // not done, gogc generation references in the object might
      // not get marked.
      // For the case of arrays, which are otherwise precisely
      // marked, we need to dirty the entire array, not just its head.
      if (is_obj_array) {
        // The [par_]mark_range() method expects mr.end() below to
        // be aligned to the granularity of a bit's representation
        // in the heap. In the case of the MUT below, that's a
        // card size.
        MemRegion mr(start,
                     align_up(start + obj_size,
                              CardTable::card_size /* bytes */));
        if (par) {
          _modUnionTable.par_mark_range(mr);
        } else {
          _modUnionTable.mark_range(mr);
        }
      } else {  // not an obj array; we can just mark the head
        if (par) {
          _modUnionTable.par_mark(start);
        } else {
          _modUnionTable.mark(start);
        }
      }
    }
  }
}

oop goConcurrentMarkSweepGeneration::promote(oop obj, size_t obj_size) {
  assert(obj_size == (size_t)obj->size(), "bad obj_size passed in");
  // allocate, copy and if necessary update promoinfo --
  // delegate to underlying space.
  assert_lock_strong(freelistLock());

#ifndef PRODUCT
  if (gogcHeap::heap()->promotion_should_fail()) {
    return NULL;
  }
#endif  // #ifndef PRODUCT

  oop res = _gogcSpace->promote(obj, obj_size);
  if (res == NULL) {
    // expand and retry
    size_t s = _gogcSpace->expansionSpaceRequired(obj_size);  // HeapWords
    expand_for_gc_cause(s*HeapWordSize, MinHeapDeltaBytes, gogcExpansionCause::_satisfy_promotion);
    // Since this is the old generation, we don't try to promote
    // into a more senior generation.
    res = _gogcSpace->promote(obj, obj_size);
  }
  if (res != NULL) {
    // See comment in allocate() about when objects should
    // be allocated live.
    assert(oopDesc::is_oop(obj), "Will dereference klass pointer below");
    collector()->promoted(false,           // Not parallel
                          (HeapWord*)res, obj->is_objArray(), obj_size);
    // promotion counters
    NOT_PRODUCT(
      _numObjectsPromoted++;
      _numWordsPromoted +=
        (int)(CompactibleFreeListSpace::adjustObjectSize(obj->size()));
    )
  }
  return res;
}

goConcurrentMarkSweepGeneration::par_promote(int thread_num,
                                           oop old, markOop m,
                                           size_t word_sz) {
#ifndef PRODUCT
  if (gogcHeap::heap()->promotion_should_fail()) {
    return NULL;
  }
#endif  // #ifndef PRODUCT

  gogcParGCThreadState* ps = _par_gc_thread_states[thread_num];
  PromotionInfo* promoInfo = &ps->promo;
   if (promoInfo->tracking() && !promoInfo->ensure_spooling_space()) {
    if (!expand_and_ensure_spooling_space(promoInfo)) {
      return NULL;
    }
  }
  assert(!promoInfo->tracking() || promoInfo->has_spooling_space(), "Control point invariant");
  const size_t alloc_sz = CompactibleFreeListSpace::adjustObjectSize(word_sz);
  HeapWord* obj_ptr = ps->lab.alloc(alloc_sz);
  if (obj_ptr == NULL) {
     obj_ptr = expand_and_par_lab_allocate(ps, alloc_sz);
     if (obj_ptr == NULL) {
       return NULL;
     }
  }
  oop obj = oop(obj_ptr);
  OrderAccess::storestore();
  assert(obj->klass_or_null() == NULL, "Object should be uninitialized here.");
  assert(!((FreeChunk*)obj_ptr)->is_free(), "Error, block will look free but show wrong size");
  HeapWord* old_ptr = (HeapWord*)old;
  obj->set_mark_raw(m);
  assert(obj->klass_or_null() == NULL, "Object should be uninitialized here.");
  assert(!((FreeChunk*)obj_ptr)->is_free(), "Error, block will look free but show wrong size");
  OrderAccess::storestore();

  if (UseCompressedClassPointers) {
    obj->set_klass_gap(old->klass_gap());
  }
  if (word_sz > (size_t)oopDesc::header_size()) {
    Copy::aligned_disjoint_words(old_ptr + oopDesc::header_size(),
                                 obj_ptr + oopDesc::header_size(),
                                 word_sz - oopDesc::header_size());
  }

    if (promoInfo->tracking()) {
    promoInfo->track((PromotedObject*)obj, old->klass());
  }
  assert(obj->klass_or_null() == NULL, "Object should be uninitialized here.");
  assert(!((FreeChunk*)obj_ptr)->is_free(), "Error, block will look free but show wrong size");
  assert(oopDesc::is_oop(old), "Will use and dereference old klass ptr below");

  
  OrderAccess::storestore();
  obj->set_klass(old->klass());
  assert(oopDesc::is_oop(obj) && obj->size() == (int)word_sz, "Error, incorrect size computed for promoted object");

  collector()->promoted(true,          // parallel
                        obj_ptr, old->is_objArray(), word_sz);

  NOT_PRODUCT(
    Atomic::inc(&_numObjectsPromoted);
    Atomic::add(alloc_sz, &_numWordsPromoted);
  )

  return obj;
}

void
goConcurrentMarkSweepGeneration::
par_promote_alloc_done(int thread_num) {
  gogcParGCThreadState* ps = _par_gc_thread_states[thread_num];
  ps->lab.retire(thread_num);
}

void
goConcurrentMarkSweepGeneration::
par_oop_since_save_marks_iterate_done(int thread_num) {
  gogcParGCThreadState* ps = _par_gc_thread_states[thread_num];
  ParScanWithoutBarrierClosure* dummy_cl = NULL;
  ps->promo.promoted_oops_iterate(dummy_cl);

  ps->promo.stopTrackingPromotions();
}

bool goConcurrentMarkSweepGeneration::should_collect(bool   full,
                                                   size_t size,
                                                   bool   tlab)
{
  return full || should_allocate(size, tlab); // FIX ME !!!
 
}

bool gogcCollector::shouldConcurrentCollect() {
  LogTarget(Trace, gc) log;

  if (_full_gc_requested) {
    log.print("gogcCollector: collect because of explicit  gc request (or GCLocker)");
    return true;
  }

  FreelistLocker x(this);
  // ------------------------------------------------------------------
  // Print out lots of information which affects the initiation of
  // a collection.
  if (log.is_enabled() && stats().valid()) {
    log.print("gogcCollector shouldConcurrentCollect: ");

    LogStream out(log);
    stats().print_on(&out);

    log.print("time_until_gogc_gen_full %3.7f", stats().time_until_gogc_gen_full());
    log.print("free=" SIZE_FORMAT, _gogcGen->free());
    log.print("contiguous_available=" SIZE_FORMAT, _gogcGen->contiguous_available());
    log.print("promotion_rate=%g", stats().promotion_rate());
    log.print("gogc_allocation_rate=%g", stats().gogc_allocation_rate());
    log.print("occupancy=%3.7f", _gogcGen->occupancy());
    log.print("initiatingOccupancy=%3.7f", _gogcGen->initiating_occupancy());
    log.print("gogc_time_since_begin=%3.7f", stats().gogc_time_since_begin());
    log.print("gogc_time_since_end=%3.7f", stats().gogc_time_since_end());
    log.print("metadata initialized %d", MetaspaceGC::should_concurrent_collect());
  }
  // ------------------------------------------------------------------

  // If the estimated time to complete a gogc collection (gogc_duration())
  // is less than the estimated time remaining until the gogc generation
  // is full, start a collection.
  if (!UsegogcInitiatingOccupancyOnly) {
    if (stats().valid()) {
      if (stats().time_until_gogc_start() == 0.0) {
        return true;
      }
    } else {
      // We want to conservatively collect somewhat early in order
      // to try and "bootstrap" our gogc/promotion statistics;
      // this branch will not fire after the first successful gogc
      // collection because the stats should then be valid.
      if (_gogcGen->occupancy() >= _bootstrap_occupancy) {
        log.print(" gogcCollector: collect for bootstrapping statistics: occupancy = %f, boot occupancy = %f",
                  _gogcGen->occupancy(), _bootstrap_occupancy);
        return true;
      }
    }
  }

  // Otherwise, we start a collection cycle if
  // old gen want a collection cycle started. Each may use
  // an appropriate criterion for making this decision.
  // XXX We need to make sure that the gen expansion
  // criterion dovetails well with this. XXX NEED TO FIX THIS
  if (_gogcGen->should_concurrent_collect()) {
    log.print("gogc old gen initiated");
    return true;
  }

  // We start a collection if we believe an incremental collection may fail;
  // this is not likely to be productive in practice because it's probably too
  // late anyway.
  gogcHeap* heap = gogcHeap::heap();
  if (heap->incremental_collection_will_fail(true /* consult_young */)) {
    log.print("gogcCollector: collect because incremental collection will fail ");
    return true;
  }

  if (MetaspaceGC::should_concurrent_collect()) {
    log.print("gogcCollector: collect for metadata allocation ");
    return true;
  }

  // gogcTriggerInterval starts a gogc cycle if enough time has passed.
  if (gogcTriggerInterval >= 0) {
    if (gogcTriggerInterval == 0) {
      // Trigger always
      return true;
    }

    // Check the gogc time since begin (we do not check the stats validity
    // as we want to be able to trigger the first gogc cycle as well)
    if (stats().gogc_time_since_begin() >= (gogcTriggerInterval / ((double) MILLIUNITS))) {
      if (stats().valid()) {
        log.print("gogcCollector: collect because of trigger interval (time since last begin %3.7f secs)",
                  stats().gogc_time_since_begin());
      } else {
        log.print("gogcCollector: collect because of trigger interval (first collection)");
      }
      return true;
    }
  }

  return false;
}

void gogcCollector::set_did_compact(bool v) { _gogcGen->set_did_compact(v); }

// Clear _expansion_cause fields of constituent generations
void gogcCollector::clear_expansion_cause() {
  _gogcGen->clear_expansion_cause();
}

// We should be conservative in starting a collection cycle.  To
// start too eagerly runs the risk of collecting too often in the
// extreme.  To collect too rarely falls back on full collections,
// which works, even if not optimum in terms of concurrent work.
// As a work around for too eagerly collecting, use the flag
// UsegogcInitiatingOccupancyOnly.  This also has the advantage of
// giving the user an easily understandable way of controlling the
// collections.
// We want to start a new collection cycle if any of the following
// conditions hold:
// . our current occupancy exceeds the configured initiating occupancy
//   for this generation, or
// . we recently needed to expand this space and have not, since that
//   expansion, done a collection of this generation, or
// . the underlying space believes that it may be a good idea to initiate
//   a concurrent collection (this may be based on criteria such as the
//   following: the space uses linear allocation and linear allocation is
//   going to fail, or there is believed to be excessive fragmentation in
//   the generation, etc... or ...
// [.(currently done by gogcCollector::shouldConcurrentCollect() only for
//   the case of the old generation; see CR 6543076):
//   we may be approaching a point at which allocation requests may fail because
//   we will be out of sufficient free space given allocation rate estimates.]
bool goConcurrentMarkSweepGeneration::should_concurrent_collect() const {

  assert_lock_strong(freelistLock());
  if (occupancy() > initiating_occupancy()) {
    log_trace(gc)(" %s: collect because of occupancy %f / %f  ",
                  short_name(), occupancy(), initiating_occupancy());
    return true;
  }
  if (UsegogcInitiatingOccupancyOnly) {
    return false;
  }
  if (expansion_cause() == gogcExpansionCause::_satisfy_allocation) {
    log_trace(gc)(" %s: collect because expanded for allocation ", short_name());
    return true;
  }
  return false;
}

void goConcurrentMarkSweepGeneration::collect(bool   full,
                                            bool   clear_all_soft_refs,
                                            size_t size,
                                            bool   tlab)
{
  collector()->collect(full, clear_all_soft_refs, size, tlab);
}

void gogcCollector::collect(bool   full,
                           bool   clear_all_soft_refs,
                           size_t size,
                           bool   tlab)
{
  // The following "if" branch is present for defensive reasons.
  // In the current uses of this interface, it can be replaced with:
  // assert(!GCLocker.is_active(), "Can't be called otherwise");
  // But I am not placing that assert here to allow future
  // generality in invoking this interface.
  if (GCLocker::is_active()) {
    // A consistency test for GCLocker
    assert(GCLocker::needs_gc(), "Should have been set already");
    // Skip this foreground collection, instead
    // expanding the heap if necessary.
    // Need the free list locks for the call to free() in compute_new_size()
    compute_new_size();
    return;
  }
  acquire_control_and_collect(full, clear_all_soft_refs);
}

void gogcCollector::request_full_gc(unsigned int full_gc_count, GCCause::Cause cause) {
  gogcHeap* heap = gogcHeap::heap();
  unsigned int gc_count = heap->total_full_collections();
  if (gc_count == full_gc_count) {
    MutexLocker y(CGC_lock, Mutex::_no_safepoint_check_flag);
    _full_gc_requested = true;
    _full_gc_cause = cause;
    CGC_lock->notify();   // nudge gogc thread
  } else {
    assert(gc_count > full_gc_count, "Error: causal loop");
  }
}

bool gogcCollector::is_external_interruption() {
  GCCause::Cause cause = gogcHeap::heap()->gc_cause();
  return GCCause::is_user_requested_gc(cause) ||
         GCCause::is_serviceability_requested_gc(cause);
}

void gogcCollector::report_concurrent_mode_interruption() {
  if (is_external_interruption()) {
    log_debug(gc)("Concurrent mode interrupted");
  } else {
    log_debug(gc)("Concurrent mode failure");
    _gc_tracer_cm->report_concurrent_mode_failure();
  }
}


// The foreground and background collectors need to coordinate in order
// to make sure that they do not mutually interfere with gogc collections.
// When a background collection is active,
// the foreground collector may need to take over (preempt) and
// synchronously complete an ongoing collection. Depending on the
// frequency of the background collections and the heap usage
// of the application, this preemption can be seldom or frequent.
// There are only certain
// points in the background collection that the "collection-baton"
// can be passed to the foreground collector.
//
// The foreground collector will wait for the baton before
// starting any part of the collection.  The foreground collector
// will only wait at one location.
//
// The background collector will yield the baton before starting a new
// phase of the collection (e.g., before initial marking, marking from roots,
// precleaning, final re-mark, sweep etc.)  This is normally done at the head
// of the loop which switches the phases. The background collector does some
// of the phases (initial mark, final re-mark) with the world stopped.
// Because of locking involved in stopping the world,
// the foreground collector should not block waiting for the background
// collector when it is doing a stop-the-world phase.  The background
// collector will yield the baton at an additional point just before
// it enters a stop-the-world phase.  Once the world is stopped, the
// background collector checks the phase of the collection.  If the
// phase has not changed, it proceeds with the collection.  If the
// phase has changed, it skips that phase of the collection.  See
// the comments on the use of the Heap_lock in collect_in_background().
//
// Variable used in baton passing.
//   _foregroundGCIsActive - Set to true by the foreground collector when
//      it wants the baton.  The foreground clears it when it has finished
//      the collection.
//   _foregroundGCShouldWait - Set to true by the background collector
//        when it is running.  The foreground collector waits while
//      _foregroundGCShouldWait is true.
//  CGC_lock - monitor used to protect access to the above variables
//      and to notify the foreground and background collectors.
//  _collectorState - current state of the gogc collection.
//
// The foreground collector
//   acquires the CGC_lock
//   sets _foregroundGCIsActive
//   waits on the CGC_lock for _foregroundGCShouldWait to be false
//     various locks acquired in preparation for the collection
//     are released so as not to block the background collector
//     that is in the midst of a collection
//   proceeds with the collection
//   clears _foregroundGCIsActive
//   returns
//
// The background collector in a loop iterating on the phases of the
//      collection
//   acquires the CGC_lock
//   sets _foregroundGCShouldWait
//   if _foregroundGCIsActive is set
//     clears _foregroundGCShouldWait, notifies _CGC_lock
//     waits on _CGC_lock for _foregroundGCIsActive to become false
//     and exits the loop.
//   otherwise
//     proceed with that phase of the collection
//     if the phase is a stop-the-world phase,
//       yield the baton once more just before enqueueing
//       the stop-world gogc operation (executed by the VM thread).
//   returns after all phases of the collection are done
//

void gogcCollector::acquire_control_and_collect(bool full,
        bool clear_all_soft_refs) {
  assert(SafepointSynchronize::is_at_safepoint(), "should be at safepoint");
  assert(!Thread::current()->is_ConcurrentGC_thread(),
         "shouldn't try to acquire control from self!");

  // Start the protocol for acquiring control of the
  // collection from the background collector (aka gogc thread).
  assert(ConcurrentMarkSweepThread::vm_thread_has_gogc_token(),
         "VM thread should have gogc token");
  // Remember the possibly interrupted state of an ongoing
  // concurrent collection
  CollectorState first_state = _collectorState;

  // Signal to a possibly ongoing concurrent collection that
  // we want to do a foreground collection.
  _foregroundGCIsActive = true;

  // release locks and wait for a notify from the background collector
  // releasing the locks in only necessary for phases which
  // do yields to improve the granularity of the collection.
  assert_lock_strong(bitMapLock());
  // We need to lock the Free list lock for the space that we are
  // currently collecting.
  assert(haveFreelistLocks(), "Must be holding free list locks");
  bitMapLock()->unlock();
  releaseFreelistLocks();
  {
    MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
    if (_foregroundGCShouldWait) {
      // We are going to be waiting for action for the gogc thread;
      // it had better not be gone (for instance at shutdown)!
      assert(ConcurrentMarkSweepThread::gogct() != NULL && !ConcurrentMarkSweepThread::gogct()->has_terminated(),
             "gogc thread must be running");
      // Wait here until the background collector gives us the go-ahead
      ConcurrentMarkSweepThread::clear_gogc_flag(
        ConcurrentMarkSweepThread::gogc_vm_has_token);  // release token
      // Get a possibly blocked gogc thread going:
      //   Note that we set _foregroundGCIsActive true above,
      //   without protection of the CGC_lock.
      CGC_lock->notify();
      assert(!ConcurrentMarkSweepThread::vm_thread_wants_gogc_token(),
             "Possible deadlock");
      while (_foregroundGCShouldWait) {
        // wait for notification
        CGC_lock->wait_without_safepoint_check();
        // Possibility of delay/starvation here, since gogc token does
        // not know to give priority to VM thread? Actually, i think
        // there wouldn't be any delay/starvation, but the proof of
        // that "fact" (?) appears non-trivial. XXX 20011219YSR
      }
      ConcurrentMarkSweepThread::set_gogc_flag(
        ConcurrentMarkSweepThread::gogc_vm_has_token);
    }
  }
  // The gogc_token is already held.  Get back the other locks.
  assert(ConcurrentMarkSweepThread::vm_thread_has_gogc_token(),
         "VM thread should have gogc token");
  getFreelistLocks();
  bitMapLock()->lock_without_safepoint_check();
  log_debug(gc, state)("gogc foreground collector has asked for control " INTPTR_FORMAT " with first state %d",
                       p2i(Thread::current()), first_state);
  log_debug(gc, state)("    gets control with state %d", _collectorState);

  // Inform gogc gen if this was due to partial collection failing.
  // The gogc gen may use this fact to determine its expansion policy.
  gogcHeap* heap = gogcHeap::heap();
  if (heap->incremental_collection_will_fail(false /* don't consult_young */)) {
    assert(!_gogcGen->incremental_collection_failed(),
           "Should have been noticed, reacted to and cleared");
    _gogcGen->set_incremental_collection_failed();
  }

  if (first_state > Idling) {
    report_concurrent_mode_interruption();
  }

  set_did_compact(true);

  // If the collection is being acquired from the background
  // collector, there may be references on the discovered
  // references lists.  Abandon those references, since some
  // of them may have become unreachable after concurrent
  // discovery; the STW compacting collector will redo discovery
  // more precisely, without being subject to floating garbage.
  // Leaving otherwise unreachable references in the discovered
  // lists would require special handling.
  ref_processor()->disable_discovery();
  ref_processor()->abandon_partial_discovery();
  ref_processor()->verify_no_references_recorded();

  if (first_state > Idling) {
    save_heap_summary();
  }

  do_compaction_work(clear_all_soft_refs);

  // Has the GC time limit been exceeded?
  size_t max_eden_size = _young_gen->max_eden_size();
  GCCause::Cause gc_cause = heap->gc_cause();
  size_policy()->check_gc_overhead_limit(_young_gen->eden()->used(),
                                         _gogcGen->max_capacity(),
                                         max_eden_size,
                                         full,
                                         gc_cause,
                                         heap->soft_ref_policy());

  // Reset the expansion cause, now that we just completed
  // a collection cycle.
  clear_expansion_cause();
  _foregroundGCIsActive = false;
  return;
}

// Resize the tenured generation
// after obtaining the free list locks for the
// two generations.
void gogcCollector::compute_new_size() {
  assert_locked_or_safepoint(Heap_lock);
  FreelistLocker z(this);
  MetaspaceGC::compute_new_size();
  _gogcGen->compute_new_size_free_list();
}

// A work method used by the foreground collector to do
// a mark-sweep-compact.
void gogcCollector::do_compaction_work(bool clear_all_soft_refs) {
  gogcHeap* heap = gogcHeap::heap();

  STWGCTimer* gc_timer = GenMarkSweep::gc_timer();
  gc_timer->register_gc_start();

  SerialOldTracer* gc_tracer = GenMarkSweep::gc_tracer();
  gc_tracer->report_gc_start(heap->gc_cause(), gc_timer->gc_start());

  heap->pre_full_gc_dump(gc_timer);

  GCTraceTime(Trace, gc, phases) t("gogc:MSC");

  // Temporarily widen the span of the weak reference processing to
  // the entire heap.
  MemRegion new_span(gogcHeap::heap()->reserved_region());
  ReferenceProcessorSpanMutator rp_mut_span(ref_processor(), new_span);
  // Temporarily, clear the "is_alive_non_header" field of the
  // reference processor.
  ReferenceProcessorIsAliveMutator rp_mut_closure(ref_processor(), NULL);
  // Temporarily make reference _processing_ single threaded (non-MT).
  ReferenceProcessorMTProcMutator rp_mut_mt_processing(ref_processor(), false);
  // Temporarily make refs discovery atomic
  ReferenceProcessorAtomicMutator rp_mut_atomic(ref_processor(), true);
  // Temporarily make reference _discovery_ single threaded (non-MT)
  ReferenceProcessorMTDiscoveryMutator rp_mut_discovery(ref_processor(), false);

  ref_processor()->set_enqueuing_is_done(false);
  ref_processor()->enable_discovery();
  ref_processor()->setup_policy(clear_all_soft_refs);
  // If an asynchronous collection finishes, the _modUnionTable is
  // all clear.  If we are assuming the collection from an asynchronous
  // collection, clear the _modUnionTable.
  assert(_collectorState != Idling || _modUnionTable.isAllClear(),
    "_modUnionTable should be clear if the baton was not passed");
  _modUnionTable.clear_all();
  assert(_collectorState != Idling || _ct->cld_rem_set()->mod_union_is_clear(),
    "mod union for klasses should be clear if the baton was passed");
  _ct->cld_rem_set()->clear_mod_union();


  // We must adjust the allocation statistics being maintained
  // in the free list space. We do so by reading and clearing
  // the sweep timer and updating the block flux rate estimates below.
  assert(!_intra_sweep_timer.is_active(), "_intra_sweep_timer should be inactive");
  if (_inter_sweep_timer.is_active()) {
    _inter_sweep_timer.stop();
    // Note that we do not use this sample to update the _inter_sweep_estimate.
    _gogcGen->gogcSpace()->beginSweepFLCensus((float)(_inter_sweep_timer.seconds()),
                                            _inter_sweep_estimate.padded_average(),
                                            _intra_sweep_estimate.padded_average());
  }

  GenMarkSweep::invoke_at_safepoint(ref_processor(), clear_all_soft_refs);
  #ifdef ASSERT
    CompactibleFreeListSpace* gogc_space = _gogcGen->gogcSpace();
    size_t free_size = gogc_space->free();
    assert(free_size ==
           pointer_delta(gogc_space->end(), gogc_space->compaction_top())
           * HeapWordSize,
      "All the free space should be compacted into one chunk at top");
    assert(gogc_space->dictionary()->total_chunk_size(
                                      debug_only(gogc_space->freelistLock())) == 0 ||
           gogc_space->totalSizeInIndexedFreeLists() == 0,
      "All the free space should be in a single chunk");
    size_t num = gogc_space->totalCount();
    assert((free_size == 0 && num == 0) ||
           (free_size > 0  && (num == 1 || num == 2)),
         "There should be at most 2 free chunks after compaction");
  #endif // ASSERT
  _collectorState = Resetting;
  assert(_restart_addr == NULL,
         "Should have been NULL'd before baton was passed");
  reset_stw();
  _gogcGen->reset_after_compaction();
  _concurrent_cycles_since_last_unload = 0;

  // Clear any data recorded in the PLAB chunk arrays.
  if (_survivor_plab_array != NULL) {
    reset_survivor_plab_arrays();
  }

  // Adjust the per-size allocation stats for the next epoch.
  _gogcGen->gogcSpace()->endSweepFLCensus(sweep_count() /* fake */);
  // Restart the "inter sweep timer" for the next epoch.
  _inter_sweep_timer.reset();
  _inter_sweep_timer.start();

  // No longer a need to do a concurrent collection for Metaspace.
  MetaspaceGC::set_should_concurrent_collect(false);

  heap->post_full_gc_dump(gc_timer);

  gc_timer->register_gc_end();

  gc_tracer->report_gc_end(gc_timer->gc_end(), gc_timer->time_partitions());

  // For a mark-sweep-compact, compute_new_size() will be called
  // in the heap's do_collection() method.
}

void gogcCollector::print_eden_and_survivor_chunk_arrays() {
  Log(gc, heap) log;
  if (!log.is_trace()) {
    return;
  }

  ContiguousSpace* eden_space = _young_gen->eden();
  ContiguousSpace* from_space = _young_gen->from();
  ContiguousSpace* to_space   = _young_gen->to();
  // Eden
  if (_eden_chunk_array != NULL) {
    log.trace("eden " PTR_FORMAT "-" PTR_FORMAT "-" PTR_FORMAT "(" SIZE_FORMAT ")",
              p2i(eden_space->bottom()), p2i(eden_space->top()),
              p2i(eden_space->end()), eden_space->capacity());
    log.trace("_eden_chunk_index=" SIZE_FORMAT ", _eden_chunk_capacity=" SIZE_FORMAT,
              _eden_chunk_index, _eden_chunk_capacity);
    for (size_t i = 0; i < _eden_chunk_index; i++) {
      log.trace("_eden_chunk_array[" SIZE_FORMAT "]=" PTR_FORMAT, i, p2i(_eden_chunk_array[i]));
    }
  }
  // Survivor
  if (_survivor_chunk_array != NULL) {
    log.trace("survivor " PTR_FORMAT "-" PTR_FORMAT "-" PTR_FORMAT "(" SIZE_FORMAT ")",
              p2i(from_space->bottom()), p2i(from_space->top()),
              p2i(from_space->end()), from_space->capacity());
    log.trace("_survivor_chunk_index=" SIZE_FORMAT ", _survivor_chunk_capacity=" SIZE_FORMAT,
              _survivor_chunk_index, _survivor_chunk_capacity);
    for (size_t i = 0; i < _survivor_chunk_index; i++) {
      log.trace("_survivor_chunk_array[" SIZE_FORMAT "]=" PTR_FORMAT, i, p2i(_survivor_chunk_array[i]));
    }
  }
}

void gogcCollector::getFreelistLocks() const {
  // Get locks for all free lists in all generations that this
  // collector is responsible for
  _gogcGen->freelistLock()->lock_without_safepoint_check();
}

void gogcCollector::releaseFreelistLocks() const {
  // Release locks for all free lists in all generations that this
  // collector is responsible for
  _gogcGen->freelistLock()->unlock();
}

bool gogcCollector::haveFreelistLocks() const {
  // Check locks for all free lists in all generations that this
  // collector is responsible for
  assert_lock_strong(_gogcGen->freelistLock());
  PRODUCT_ONLY(ShouldNotReachHere());
  return true;
}

// A utility class that is used by the gogc collector to
// temporarily "release" the foreground collector from its
// usual obligation to wait for the background collector to
// complete an ongoing phase before proceeding.
class ReleaseForegroundGC: public StackObj {
 private:
  gogcCollector* _c;
 public:
  ReleaseForegroundGC(gogcCollector* c) : _c(c) {
    assert(_c->_foregroundGCShouldWait, "Else should not need to call");
    MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
    // allow a potentially blocked foreground collector to proceed
    _c->_foregroundGCShouldWait = false;
    if (_c->_foregroundGCIsActive) {
      CGC_lock->notify();
    }
    assert(!ConcurrentMarkSweepThread::gogc_thread_has_gogc_token(),
           "Possible deadlock");
  }

  ~ReleaseForegroundGC() {
    assert(!_c->_foregroundGCShouldWait, "Usage protocol violation?");
    MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
    _c->_foregroundGCShouldWait = true;
  }
};

void gogcCollector::collect_in_background(GCCause::Cause cause) {
  assert(Thread::current()->is_ConcurrentGC_thread(),
    "A gogc asynchronous collection is only allowed on a gogc thread.");

  gogcHeap* heap = gogcHeap::heap();
  {
    MutexLocker hl(Heap_lock, Mutex::_no_safepoint_check_flag);
    FreelistLocker fll(this);
    MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
    if (_foregroundGCIsActive) {
      // The foreground collector is. Skip this
      // background collection.
      assert(!_foregroundGCShouldWait, "Should be clear");
      return;
    } else {
      assert(_collectorState == Idling, "Should be idling before start.");
      _collectorState = InitialMarking;
      register_gc_start(cause);
      // Reset the expansion cause, now that we are about to begin
      // a new cycle.
      clear_expansion_cause();

      // Clear the MetaspaceGC flag since a concurrent collection
      // is starting but also clear it after the collection.
      MetaspaceGC::set_should_concurrent_collect(false);
    }
    // Decide if we want to enable class unloading as part of the
    // ensuing concurrent GC cycle.
    update_should_unload_classes();
    _full_gc_requested = false;           // acks all outstanding full gc requests
    _full_gc_cause = GCCause::_no_gc;
    // Signal that we are about to start a collection
    heap->increment_total_full_collections();  // ... starting a collection cycle
    _collection_count_start = heap->total_full_collections();
  }

  size_t prev_used = _gogcGen->used();

  // The change of the collection state is normally done at this level;
  // the exceptions are phases that are executed while the world is
  // stopped.  For those phases the change of state is done while the
  // world is stopped.  For baton passing purposes this allows the
  // background collector to finish the phase and change state atomically.
  // The foreground collector cannot wait on a phase that is done
  // while the world is stopped because the foreground collector already
  // has the world stopped and would deadlock.
  while (_collectorState != Idling) {
    log_debug(gc, state)("Thread " INTPTR_FORMAT " in gogc state %d",
                         p2i(Thread::current()), _collectorState);
    // The foreground collector
    //   holds the Heap_lock throughout its collection.
    //   holds the gogc token (but not the lock)
    //     except while it is waiting for the background collector to yield.
    //
    // The foreground collector should be blocked (not for long)
    //   if the background collector is about to start a phase
    //   executed with world stopped.  If the background
    //   collector has already started such a phase, the
    //   foreground collector is blocked waiting for the
    //   Heap_lock.  The stop-world phases (InitialMarking and FinalMarking)
    //   are executed in the VM thread.
    //
    // The locking order is
    //   PendingListLock (PLL)  -- if applicable (FinalMarking)
    //   Heap_lock  (both this & PLL locked in VM_gogc_Operation::prologue())
    //   gogc token  (claimed in
    //                stop_world_and_do() -->
    //                  safepoint_synchronize() -->
    //                    gogcThread::synchronize())

    {
      // Check if the FG collector wants us to yield.
      gogcTokenSync x(true); // is gogc thread
      if (waitForForegroundGC()) {
        // We yielded to a foreground GC, nothing more to be
        // done this round.
        assert(_foregroundGCShouldWait == false, "We set it to false in "
               "waitForForegroundGC()");
        log_debug(gc, state)("gogc Thread " INTPTR_FORMAT " exiting collection gogc state %d",
                             p2i(Thread::current()), _collectorState);
        return;
      } else {
        // The background collector can run but check to see if the
        // foreground collector has done a collection while the
        // background collector was waiting to get the CGC_lock
        // above.  If yes, break so that _foregroundGCShouldWait
        // is cleared before returning.
        if (_collectorState == Idling) {
          break;
        }
      }
    }

    assert(_foregroundGCShouldWait, "Foreground collector, if active, "
      "should be waiting");

    switch (_collectorState) {
      case InitialMarking:
        {
          ReleaseForegroundGC x(this);
          stats().record_gogc_begin();
          VM_gogc_Initial_Mark initial_mark_op(this);
          VMThread::execute(&initial_mark_op);
        }
        // The collector state may be any legal state at this point
        // since the background collector may have yielded to the
        // foreground collector.
        break;
      case Marking:
        // initial marking in checkpointRootsInitialWork has been completed
        if (markFromRoots()) { // we were successful
          assert(_collectorState == Precleaning, "Collector state should "
            "have changed");
        } else {
          assert(_foregroundGCIsActive, "Internal state inconsistency");
        }
        break;
      case Precleaning:
        // marking from roots in markFromRoots has been completed
        preclean();
        assert(_collectorState == AbortablePreclean ||
               _collectorState == FinalMarking,
               "Collector state should have changed");
        break;
      case AbortablePreclean:
        abortable_preclean();
        assert(_collectorState == FinalMarking, "Collector state should "
          "have changed");
        break;
      case FinalMarking:
        {
          ReleaseForegroundGC x(this);

          VM_gogc_Final_Remark final_remark_op(this);
          VMThread::execute(&final_remark_op);
        }
        assert(_foregroundGCShouldWait, "block post-condition");
        break;
      case Sweeping:
        // final marking in checkpointRootsFinal has been completed
        sweep();
        assert(_collectorState == Resizing, "Collector state change "
          "to Resizing must be done under the free_list_lock");

      case Resizing: {
        // Sweeping has been completed...
        // At this point the background collection has completed.
        // Don't move the call to compute_new_size() down
        // into code that might be executed if the background
        // collection was preempted.
        {
          ReleaseForegroundGC x(this);   // unblock FG collection
          MutexLocker         y(Heap_lock, Mutex::_no_safepoint_check_flag);
          gogcTokenSync        z(true);   // not strictly needed.
          if (_collectorState == Resizing) {
            compute_new_size();
            save_heap_summary();
            _collectorState = Resetting;
          } else {
            assert(_collectorState == Idling, "The state should only change"
                   " because the foreground collector has finished the collection");
          }
        }
        break;
      }
      case Resetting:
        // gogc heap resizing has been completed
        reset_concurrent();
        assert(_collectorState == Idling, "Collector state should "
          "have changed");

        MetaspaceGC::set_should_concurrent_collect(false);

        stats().record_gogc_end();
        // Don't move the concurrent_phases_end() and compute_new_size()
        // calls to here because a preempted background collection
        // has it's state set to "Resetting".
        break;
      case Idling:
      default:
        ShouldNotReachHere();
        break;
    }
    log_debug(gc, state)("  Thread " INTPTR_FORMAT " done - next gogc state %d",
                         p2i(Thread::current()), _collectorState);
    assert(_foregroundGCShouldWait, "block post-condition");
  }

  // Should this be in gc_epilogue?
  heap->counters()->update_counters();

  {
    // Clear _foregroundGCShouldWait and, in the event that the
    // foreground collector is waiting, notify it, before
    // returning.
    MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
    _foregroundGCShouldWait = false;
    if (_foregroundGCIsActive) {
      CGC_lock->notify();
    }
    assert(!ConcurrentMarkSweepThread::gogc_thread_has_gogc_token(),
           "Possible deadlock");
  }
  log_debug(gc, state)("gogc Thread " INTPTR_FORMAT " exiting collection gogc state %d",
                       p2i(Thread::current()), _collectorState);
  log_info(gc, heap)("Old: " SIZE_FORMAT "K->" SIZE_FORMAT "K("  SIZE_FORMAT "K)",
                     prev_used / K, _gogcGen->used()/K, _gogcGen->capacity() /K);
}

void gogcCollector::register_gc_start(GCCause::Cause cause) {
  _gogc_start_registered = true;
  _gc_timer_cm->register_gc_start();
  _gc_tracer_cm->report_gc_start(cause, _gc_timer_cm->gc_start());
}

void gogcCollector::register_gc_end() {
  if (_gogc_start_registered) {
    report_heap_summary(GCWhen::AfterGC);

    _gc_timer_cm->register_gc_end();
    _gc_tracer_cm->report_gc_end(_gc_timer_cm->gc_end(), _gc_timer_cm->time_partitions());
    _gogc_start_registered = false;
  }
}

void gogcCollector::save_heap_summary() {
  gogcHeap* heap = gogcHeap::heap();
  _last_heap_summary = heap->create_heap_summary();
  _last_metaspace_summary = heap->create_metaspace_summary();
}

void gogcCollector::report_heap_summary(GCWhen::Type when) {
  _gc_tracer_cm->report_gc_heap_summary(when, _last_heap_summary);
  _gc_tracer_cm->report_metaspace_summary(when, _last_metaspace_summary);
}

bool gogcCollector::waitForForegroundGC() {
  bool res = false;
  assert(ConcurrentMarkSweepThread::gogc_thread_has_gogc_token(),
         "gogc thread should have gogc token");
  // Block the foreground collector until the
  // background collectors decides whether to
  // yield.
  MutexLocker x(CGC_lock, Mutex::_no_safepoint_check_flag);
  _foregroundGCShouldWait = true;
  if (_foregroundGCIsActive) {
    // The background collector yields to the
    // foreground collector and returns a value
    // indicating that it has yielded.  The foreground
    // collector can proceed.
    res = true;
    _foregroundGCShouldWait = false;
    ConcurrentMarkSweepThread::clear_gogc_flag(
      ConcurrentMarkSweepThread::gogc_gogc_has_token);
    ConcurrentMarkSweepThread::set_gogc_flag(
      ConcurrentMarkSweepThread::gogc_gogc_wants_token);
    // Get a possibly blocked foreground thread going
    CGC_lock->notify();
    log_debug(gc, state)("gogc Thread " INTPTR_FORMAT " waiting at gogc state %d",
                         p2i(Thread::current()), _collectorState);
    while (_foregroundGCIsActive) {
      CGC_lock->wait_without_safepoint_check();
    }
    ConcurrentMarkSweepThread::set_gogc_flag(
      ConcurrentMarkSweepThread::gogc_gogc_has_token);
    ConcurrentMarkSweepThread::clear_gogc_flag(
      ConcurrentMarkSweepThread::gogc_gogc_wants_token);
  }
  log_debug(gc, state)("gogc Thread " INTPTR_FORMAT " continuing at gogc state %d",
                       p2i(Thread::current()), _collectorState);
  return res;
}

// Because of the need to lock the free lists and other structures in
// the collector, common to all the generations that the collector is
// collecting, we need the gc_prologues of individual gogc generations
// delegate to their collector. It may have been simpler had the
// current infrastructure allowed one to call a prologue on a
// collector. In the absence of that we have the generation's
// prologue delegate to the collector, which delegates back
// some "local" work to a worker method in the individual generations
// that it's responsible for collecting, while itself doing any
// work common to all generations it's responsible for. A similar
// comment applies to the  gc_epilogue()'s.
// The role of the variable _between_prologue_and_epilogue is to
// enforce the invocation protocol.
void gogcCollector::gc_prologue(bool full) {
  // Call gc_prologue_work() for the gogcGen
  // we are responsible for.

  // The following locking discipline assumes that we are only called
  // when the world is stopped.
  assert(SafepointSynchronize::is_at_safepoint(), "world is stopped assumption");

  // The gogcCollector prologue must call the gc_prologues for the
  // "generations" that it's responsible
  // for.

  assert(   Thread::current()->is_VM_thread()
         || (   gogcScavengeBeforeRemark
             && Thread::current()->is_ConcurrentGC_thread()),
         "Incorrect thread type for prologue execution");

  if (_between_prologue_and_epilogue) {
    // We have already been invoked; this is a gc_prologue delegation
    // from yet another gogc generation that we are responsible for, just
    // ignore it since all relevant work has already been done.
    return;
  }

  // set a bit saying prologue has been called; cleared in epilogue
  _between_prologue_and_epilogue = true;
  // Claim locks for common data structures, then call gc_prologue_work()
  // for each gogcGen.

  getFreelistLocks();   // gets free list locks on constituent spaces
  bitMapLock()->lock_without_safepoint_check();

  // Should call gc_prologue_work() for all gogc gens we are responsible for
  bool duringMarking =    _collectorState >= Marking
                         && _collectorState < Sweeping;

  // The young collections clear the modified oops state, which tells if
  // there are any modified oops in the class. The remark phase also needs
  // that information. Tell the young collection to save the union of all
  // modified klasses.
  if (duringMarking) {
    _ct->cld_rem_set()->set_accumulate_modified_oops(true);
  }

  bool registerClosure = duringMarking;

  _gogcGen->gc_prologue_work(full, registerClosure, &_modUnionClosurePar);

  if (!full) {
    stats().record_gc0_begin();
  }
}

void goConcurrentMarkSweepGeneration::gc_prologue(bool full) {

  _capacity_at_prologue = capacity();
  _used_at_prologue = used();

  // We enable promotion tracking so that card-scanning can recognize
  // which objects have been promoted during this GC and skip them.
  for (uint i = 0; i < ParallelGCThreads; i++) {
    _par_gc_thread_states[i]->promo.startTrackingPromotions();
  }

  // Delegate to gogccollector which knows how to coordinate between
  // this and any other gogc generations that it is responsible for
  // collecting.
  collector()->gc_prologue(full);
}

// This is a "private" interface for use by this generation's gogcCollector.
// Not to be called directly by any other entity (for instance,
// GenCollectedHeap, which calls the "public" gc_prologue method above).
void goConcurrentMarkSweepGeneration::gc_prologue_work(bool full,
  bool registerClosure, ModUnionClosure* modUnionClosure) {
  assert(!incremental_collection_failed(), "Shouldn't be set yet");
  assert(gogcSpace()->preconsumptionDirtyCardClosure() == NULL,
    "Should be NULL");
  if (registerClosure) {
    gogcSpace()->setPreconsumptionDirtyCardClosure(modUnionClosure);
  }
  gogcSpace()->gc_prologue();
  // Clear stat counters
  NOT_PRODUCT(
    assert(_numObjectsPromoted == 0, "check");
    assert(_numWordsPromoted   == 0, "check");
    log_develop_trace(gc, alloc)("Allocated " SIZE_FORMAT " objects, " SIZE_FORMAT " bytes concurrently",
                                 _numObjectsAllocated, _numWordsAllocated*sizeof(HeapWord));
    _numObjectsAllocated = 0;
    _numWordsAllocated   = 0;
  )
}

void gogcCollector::gc_epilogue(bool full) {
  // The following locking discipline assumes that we are only called
  // when the world is stopped.
  assert(SafepointSynchronize::is_at_safepoint(),
         "world is stopped assumption");

  // Currently the gogc epilogue (see CompactibleFreeListSpace) merely checks
  // if linear allocation blocks need to be appropriately marked to allow the
  // the blocks to be parsable. We also check here whether we need to nudge the
  // gogc collector thread to start a new cycle (if it's not already active).
  assert(   Thread::current()->is_VM_thread()
         || (   gogcScavengeBeforeRemark
             && Thread::current()->is_ConcurrentGC_thread()),
         "Incorrect thread type for epilogue execution");

  if (!_between_prologue_and_epilogue) {
    // We have already been invoked; this is a gc_epilogue delegation
    // from yet another gogc generation that we are responsible for, just
    // ignore it since all relevant work has already been done.
    return;
  }
  assert(haveFreelistLocks(), "must have freelist locks");
  assert_lock_strong(bitMapLock());

  _ct->cld_rem_set()->set_accumulate_modified_oops(false);

  _gogcGen->gc_epilogue_work(full);

  if (_collectorState == AbortablePreclean || _collectorState == Precleaning) {
    // in case sampling was not already enabled, enable it
    _start_sampling = true;
  }
  // reset _eden_chunk_array so sampling starts afresh
  _eden_chunk_index = 0;

  size_t gogc_used   = _gogcGen->gogcSpace()->used();

  // update performance counters - this uses a special version of
  // update_counters() that allows the utilization to be passed as a
  // parameter, avoiding multiple calls to used().
  //
  _gogcGen->update_counters(gogc_used);

  bitMapLock()->unlock();
  releaseFreelistLocks();

  if (!CleanChunkPoolAsync) {
    Chunk::clean_chunk_pool();
  }

  set_did_compact(false);
  _between_prologue_and_epilogue = false;  // ready for next cycle
}

void goConcurrentMarkSweepGeneration::gc_epilogue(bool full) {
  collector()->gc_epilogue(full);

  // When using ParNew, promotion tracking should have already been
  // disabled. However, the prologue (which enables promotion
  // tracking) and epilogue are called irrespective of the type of
  // GC. So they will also be called before and after Full GCs, during
  // which promotion tracking will not be explicitly disabled. So,
  // it's safer to also disable it here too (to be symmetric with
  // enabling it in the prologue).
  for (uint i = 0; i < ParallelGCThreads; i++) {
    _par_gc_thread_states[i]->promo.stopTrackingPromotions();
  }
}

void goConcurrentMarkSweepGeneration::gc_epilogue_work(bool full) {
  assert(!incremental_collection_failed(), "Should have been cleared");
  gogcSpace()->setPreconsumptionDirtyCardClosure(NULL);
  gogcSpace()->gc_epilogue();
    // Print stat counters
  NOT_PRODUCT(
    assert(_numObjectsAllocated == 0, "check");
    assert(_numWordsAllocated == 0, "check");
    log_develop_trace(gc, promotion)("Promoted " SIZE_FORMAT " objects, " SIZE_FORMAT " bytes",
                                     _numObjectsPromoted, _numWordsPromoted*sizeof(HeapWord));
    _numObjectsPromoted = 0;
    _numWordsPromoted   = 0;
  )

  // Call down the chain in contiguous_available needs the freelistLock
  // so print this out before releasing the freeListLock.
  log_develop_trace(gc)(" Contiguous available " SIZE_FORMAT " bytes ", contiguous_available());
}

#ifndef PRODUCT
bool gogcCollector::have_gogc_token() {
  Thread* thr = Thread::current();
  if (thr->is_VM_thread()) {
    return ConcurrentMarkSweepThread::vm_thread_has_gogc_token();
  } else if (thr->is_ConcurrentGC_thread()) {
    return ConcurrentMarkSweepThread::gogc_thread_has_gogc_token();
  } else if (thr->is_GC_task_thread()) {
    return ConcurrentMarkSweepThread::vm_thread_has_gogc_token() &&
           ParGCRareEvent_lock->owned_by_self();
  }
  return false;
}

// Check reachability of the given heap address in gogc generation,
// treating all other generations as roots.
bool gogcCollector::is_gogc_reachable(HeapWord* addr) {
  // We could "guarantee" below, rather than assert, but I'll
  // leave these as "asserts" so that an adventurous debugger
  // could try this in the product build provided some subset of
  // the conditions were met, provided they were interested in the
  // results and knew that the computation below wouldn't interfere
  // with other concurrent computations mutating the structures
  // being read or written.
  assert(SafepointSynchronize::is_at_safepoint(),
         "Else mutations in object graph will make answer suspect");
  assert(have_gogc_token(), "Should hold gogc token");
  assert(haveFreelistLocks(), "must hold free list locks");
  assert_lock_strong(bitMapLock());

  // Clear the marking bit map array before starting, but, just
  // for kicks, first report if the given address is already marked
  tty->print_cr("Start: Address " PTR_FORMAT " is%s marked", p2i(addr),
                _markBitMap.isMarked(addr) ? "" : " not");

  if (verify_after_remark()) {
    MutexLocker x(verification_mark_bm()->lock(), Mutex::_no_safepoint_check_flag);
    bool result = verification_mark_bm()->isMarked(addr);
    tty->print_cr("TransitiveMark: Address " PTR_FORMAT " %s marked", p2i(addr),
                  result ? "IS" : "is NOT");
    return result;
  } else {
    tty->print_cr("Could not compute result");
    return false;
  }
}
#endif

void
gogcCollector::print_on_error(outputStream* st) {
  gogcCollector* collector = goConcurrentMarkSweepGeneration::_collector;
  if (collector != NULL) {
    gogcBitMap* bitmap = &collector->_markBitMap;
    st->print_cr("Marking Bits: (gogcBitMap*) " PTR_FORMAT, p2i(bitmap));
    bitmap->print_on_error(st, " Bits: ");

    st->cr();

    gogcBitMap* mut_bitmap = &collector->_modUnionTable;
    st->print_cr("Mod Union Table: (gogcBitMap*) " PTR_FORMAT, p2i(mut_bitmap));
    mut_bitmap->print_on_error(st, " Bits: ");
  }
}

////////////////////////////////////////////////////////
// gogc Verification Support
////////////////////////////////////////////////////////
// Following the remark phase, the following invariant
// should hold -- each object in the gogc heap which is
// marked in markBitMap() should be marked in the verification_mark_bm().

class VerifyMarkedClosure: public BitMapClosure {
  gogcBitMap* _marks;
  bool       _failed;

 public:
  VerifyMarkedClosure(gogcBitMap* bm): _marks(bm), _failed(false) {}

  bool do_bit(size_t offset) {
    HeapWord* addr = _marks->offsetToHeapWord(offset);
    if (!_marks->isMarked(addr)) {
      Log(gc, verify) log;
      ResourceMark rm;
      LogStream ls(log.error());
      oop(addr)->print_on(&ls);
      log.error(" (" INTPTR_FORMAT " should have been marked)", p2i(addr));
      _failed = true;
    }
    return true;
  }

  bool failed() { return _failed; }
};

bool gogcCollector::verify_after_remark() {
  GCTraceTime(Info, gc, phases, verify) tm("Verifying gogc Marking.");
  MutexLocker ml(verification_mark_bm()->lock(), Mutex::_no_safepoint_check_flag);
  static bool init = false;

  assert(SafepointSynchronize::is_at_safepoint(),
         "Else mutations in object graph will make answer suspect");
  assert(have_gogc_token(),
         "Else there may be mutual interference in use of "
         " verification data structures");
  assert(_collectorState > Marking && _collectorState <= Sweeping,
         "Else marking info checked here may be obsolete");
  assert(haveFreelistLocks(), "must hold free list locks");
  assert_lock_strong(bitMapLock());


  // Allocate marking bit map if not already allocated
  if (!init) { // first time
    if (!verification_mark_bm()->allocate(_span)) {
      return false;
    }
    init = true;
  }

  assert(verification_mark_stack()->isEmpty(), "Should be empty");

  // Turn off refs discovery -- so we will be tracing through refs.
  // This is as intended, because by this time
  // GC must already have cleared any refs that need to be cleared,
  // and traced those that need to be marked; moreover,
  // the marking done here is not going to interfere in any
  // way with the marking information used by GC.
  NoRefDiscovery no_discovery(ref_processor());

#if COMPILER2_OR_JVMCI
  DerivedPointerTableDeactivate dpt_deact;
#endif

  // Clear any marks from a previous round
  verification_mark_bm()->clear_all();
  assert(verification_mark_stack()->isEmpty(), "markStack should be empty");
  verify_work_stacks_empty();

  gogcHeap* heap = gogcHeap::heap();
  heap->ensure_parsability(false);  // fill TLABs, but no need to retire them
  // Update the saved marks which may affect the root scans.
  heap->save_marks();

  if (gogcRemarkVerifyVariant == 1) {
    // In this first variant of verification, we complete
    // all marking, then check if the new marks-vector is
    // a subset of the gogc marks-vector.
    verify_after_remark_work_1();
  } else {
    guarantee(gogcRemarkVerifyVariant == 2, "Range checking for gogcRemarkVerifyVariant should guarantee 1 or 2");
    // In this second variant of verification, we flag an error
    // (i.e. an object reachable in the new marks-vector not reachable
    // in the gogc marks-vector) immediately, also indicating the
    // identify of an object (A) that references the unmarked object (B) --
    // presumably, a mutation to A failed to be picked up by preclean/remark?
    verify_after_remark_work_2();
  }

  return true;
}

void gogcCollector::verify_after_remark_work_1() {
  ResourceMark rm;
  HandleMark  hm;
  gogcHeap* heap = gogcHeap::heap();

  // Get a clear set of claim bits for the roots processing to work with.
  ClassLoaderDataGraph::clear_claimed_marks();

  // Mark from roots one level into gogc
  MarkRefsIntoClosure notOlder(_span, verification_mark_bm());
  heap->rem_set()->prepare_for_younger_refs_iterate(false); // Not parallel.

  {
    StrongRootsScope srs(1);

    heap->gogc_process_roots(&srs,
                           true,   // young gen as roots
                           GenCollectedHeap::ScanningOption(roots_scanning_options()),
                           should_unload_classes(),
                           &notOlder,
                           NULL);
  }

  // Now mark from the roots
  MarkFromRootsClosure markFromRootsClosure(this, _span,
    verification_mark_bm(), verification_mark_stack(),
    false /* don't yield */, true /* verifying */);
  assert(_restart_addr == NULL, "Expected pre-condition");
  verification_mark_bm()->iterate(&markFromRootsClosure);
  while (_restart_addr != NULL) {
    // Deal with stack overflow: by restarting at the indicated
    // address.
    HeapWord* ra = _restart_addr;
    markFromRootsClosure.reset(ra);
    _restart_addr = NULL;
    verification_mark_bm()->iterate(&markFromRootsClosure, ra, _span.end());
  }
  assert(verification_mark_stack()->isEmpty(), "Should have been drained");
  verify_work_stacks_empty();

  // Marking completed -- now verify that each bit marked in
  // verification_mark_bm() is also marked in markBitMap(); flag all
  // errors by printing corresponding objects.
  VerifyMarkedClosure vcl(markBitMap());
  verification_mark_bm()->iterate(&vcl);
  if (vcl.failed()) {
    Log(gc, verify) log;
    log.error("Failed marking verification after remark");
    ResourceMark rm;
    LogStream ls(log.error());
    heap->print_on(&ls);
    fatal("gogc: failed marking verification after remark");
  }
}

class VerifyCLDOopsCLDClosure : public CLDClosure {
  class VerifyCLDOopsClosure : public OopClosure {
    gogcBitMap* _bitmap;
   public:
    VerifyCLDOopsClosure(gogcBitMap* bitmap) : _bitmap(bitmap) { }
    void do_oop(oop* p)       { guarantee(*p == NULL || _bitmap->isMarked((HeapWord*) *p), "Should be marked"); }
    void do_oop(narrowOop* p) { ShouldNotReachHere(); }
  } _oop_closure;
 public:
  VerifyCLDOopsCLDClosure(gogcBitMap* bitmap) : _oop_closure(bitmap) {}
  void do_cld(ClassLoaderData* cld) {
    cld->oops_do(&_oop_closure, ClassLoaderData::_claim_none, false);
  }
};

void gogcCollector::verify_after_remark_work_2() {
  ResourceMark rm;
  HandleMark  hm;
  gogcHeap* heap = gogcHeap::heap();

  // Get a clear set of claim bits for the roots processing to work with.
  ClassLoaderDataGraph::clear_claimed_marks();

  // Mark from roots one level into gogc
  MarkRefsIntoVerifyClosure notOlder(_span, verification_mark_bm(),
                                     markBitMap());
  CLDToOopClosure cld_closure(&notOlder, ClassLoaderData::_claim_strong);

  heap->rem_set()->prepare_for_younger_refs_iterate(false); // Not parallel.

  {
    StrongRootsScope srs(1);

    heap->gogc_process_roots(&srs,
                           true,   // young gen as roots
                           GenCollectedHeap::ScanningOption(roots_scanning_options()),
                           should_unload_classes(),
                           &notOlder,
                           &cld_closure);
  }

  // Now mark from the roots
  MarkFromRootsVerifyClosure markFromRootsClosure(this, _span,
    verification_mark_bm(), markBitMap(), verification_mark_stack());
  assert(_restart_addr == NULL, "Expected pre-condition");
  verification_mark_bm()->iterate(&markFromRootsClosure);
  while (_restart_addr != NULL) {
    // Deal with stack overflow: by restarting at the indicated
    // address.
    HeapWord* ra = _restart_addr;
    markFromRootsClosure.reset(ra);
    _restart_addr = NULL;
    verification_mark_bm()->iterate(&markFromRootsClosure, ra, _span.end());
  }
  assert(verification_mark_stack()->isEmpty(), "Should have been drained");
  verify_work_stacks_empty();

  VerifyCLDOopsCLDClosure verify_cld_oops(verification_mark_bm());
  ClassLoaderDataGraph::cld_do(&verify_cld_oops);

  // Marking completed -- now verify that each bit marked in
  // verification_mark_bm() is also marked in markBitMap(); flag all
  // errors by printing corresponding objects.
  VerifyMarkedClosure vcl(markBitMap());
  verification_mark_bm()->iterate(&vcl);
  assert(!vcl.failed(), "Else verification above should not have succeeded");
}

void goConcurrentMarkSweepGeneration::save_marks() {
  // delegate to gogc space
  gogcSpace()->save_marks();
}

bool goConcurrentMarkSweepGeneration::no_allocs_since_save_marks() {
  return gogcSpace()->no_allocs_since_save_marks();
}

void
goConcurrentMarkSweepGeneration::oop_iterate(OopIterateClosure* cl) {
  if (freelistLock()->owned_by_self()) {
    Generation::oop_iterate(cl);
  } else {
    MutexLocker x(freelistLock(), Mutex::_no_safepoint_check_flag);
    Generation::oop_iterate(cl);
  }
}

void
goConcurrentMarkSweepGeneration::object_iterate(ObjectClosure* cl) {
  if (freelistLock()->owned_by_self()) {
    Generation::object_iterate(cl);
  } else {
    MutexLocker x(freelistLock(), Mutex::_no_safepoint_check_flag);
    Generation::object_iterate(cl);
  }
}

void
goConcurrentMarkSweepGeneration::safe_object_iterate(ObjectClosure* cl) {
  if (freelistLock()->owned_by_self()) {
    Generation::safe_object_iterate(cl);
  } else {
    MutexLocker x(freelistLock(), Mutex::_no_safepoint_check_flag);
    Generation::safe_object_iterate(cl);
  }
}

void
goConcurrentMarkSweepGeneration::post_compact() {
}

void
goConcurrentMarkSweepGeneration::prepare_for_verify() {
  // Fix the linear allocation blocks to look like free blocks.

  // Locks are normally acquired/released in gc_prologue/gc_epilogue, but those
  // are not called when the heap is verified during universe initialization and
  // at vm shutdown.
  if (freelistLock()->owned_by_self()) {
    gogcSpace()->prepare_for_verify();
  } else {
    MutexLocker fll(freelistLock(), Mutex::_no_safepoint_check_flag);
    gogcSpace()->prepare_for_verify();
  }
}

void
goConcurrentMarkSweepGeneration::verify() {
  // Locks are normally acquired/released in gc_prologue/gc_epilogue, but those
  // are not called when the heap is verified during universe initialization and
  // at vm shutdown.
  if (freelistLock()->owned_by_self()) {
    gogcSpace()->verify();
  } else {
    MutexLocker fll(freelistLock(), Mutex::_no_safepoint_check_flag);
    gogcSpace()->verify();
  }
}

void gogcCollector::verify() {
  _gogcGen->verify();
}

#ifndef PRODUCT
bool gogcCollector::overflow_list_is_empty() const {
  assert(_num_par_pushes >= 0, "Inconsistency");
  if (_overflow_list == NULL) {
    assert(_num_par_pushes == 0, "Inconsistency");
  }
  return _overflow_list == NULL;
}

// The methods verify_work_stacks_empty() and verify_overflow_empty()
// merely consolidate assertion checks that appear to occur together frequently.
void gogcCollector::verify_work_stacks_empty() const {
  assert(_markStack.isEmpty(), "Marking stack should be empty");
  assert(overflow_list_is_empty(), "Overflow list should be empty");
}

void gogcCollector::verify_overflow_empty() const {
  assert(overflow_list_is_empty(), "Overflow list should be empty");
  assert(no_preserved_marks(), "No preserved marks");
}
#endif // PRODUCT

// Decide if we want to enable class unloading as part of the
// ensuing concurrent GC cycle. We will collect and
// unload classes if it's the case that:
//  (a) class unloading is enabled at the command line, and
//  (b) old gen is getting really full
// NOTE: Provided there is no change in the state of the heap between
// calls to this method, it should have idempotent results. Moreover,
// its results should be monotonically increasing (i.e. going from 0 to 1,
// but not 1 to 0) between successive calls between which the heap was
// not collected. For the implementation below, it must thus rely on
// the property that concurrent_cycles_since_last_unload()
// will not decrease unless a collection cycle happened and that
// _gogcGen->is_too_full() are
// themselves also monotonic in that sense. See check_monotonicity()
// below.
void gogcCollector::update_should_unload_classes() {
  _should_unload_classes = false;
  if (gogcClassUnloadingEnabled) {
    _should_unload_classes = (concurrent_cycles_since_last_unload() >=
                              gogcClassUnloadingMaxInterval)
                           || _gogcGen->is_too_full();
  }
}

bool goConcurrentMarkSweepGeneration::is_too_full() const {
  bool res = should_concurrent_collect();
  res = res && (occupancy() > (double)gogcIsTooFullPercentage/100.0);
  return res;
}

void gogcCollector::setup_gogc_unloading_and_verification_state() {
  const  bool should_verify =   VerifyBeforeGC || VerifyAfterGC || VerifyDuringGC
                             || VerifyBeforeExit;
  const  int  rso           =   GenCollectedHeap::SO_AllCodeCache;

  // We set the proper root for this gogc cycle here.
  if (should_unload_classes()) {   // Should unload classes this cycle
    remove_root_scanning_option(rso);  // Shrink the root set appropriately
    set_verifying(should_verify);    // Set verification state for this cycle
    return;                            // Nothing else needs to be done at this time
  }

  // Not unloading classes this cycle
  assert(!should_unload_classes(), "Inconsistency!");

  // If we are not unloading classes then add SO_AllCodeCache to root
  // scanning options.
  add_root_scanning_option(rso);

  if ((!verifying() || unloaded_classes_last_cycle()) && should_verify) {
    set_verifying(true);
  } else if (verifying() && !should_verify) {
    // We were verifying, but some verification flags got disabled.
    set_verifying(false);
    // Exclude symbols, strings and code cache elements from root scanning to
    // reduce IM and RM pauses.
    remove_root_scanning_option(rso);
  }
}


#ifndef PRODUCT
HeapWord* gogcCollector::block_start(const void* p) const {
  const HeapWord* addr = (HeapWord*)p;
  if (_span.contains(p)) {
    if (_gogcGen->gogcSpace()->is_in_reserved(addr)) {
      return _gogcGen->gogcSpace()->block_start(p);
    }
  }
  return NULL;
}
#endif

HeapWord*
goConcurrentMarkSweepGeneration::expand_and_allocate(size_t word_size,
                                                   bool   tlab,
                                                   bool   parallel) {
  gogcSynchronousYieldRequest yr;
  assert(!tlab, "Can't deal with TLAB allocation");
  MutexLocker x(freelistLock(), Mutex::_no_safepoint_check_flag);
  expand_for_gc_cause(word_size*HeapWordSize, MinHeapDeltaBytes, gogcExpansionCause::_satisfy_allocation);
  if (GCExpandToAllocateDelayMillis > 0) {
    os::sleep(Thread::current(), GCExpandToAllocateDelayMillis, false);
  }
  return have_lock_and_allocate(word_size, tlab);
}

void goConcurrentMarkSweepGeneration::expand_for_gc_cause(
    size_t bytes,
    size_t expand_bytes,
    gogcExpansionCause::Cause cause)
{

  bool success = expand(bytes, expand_bytes);

  // remember why we expanded; this information is used
  // by shouldConcurrentCollect() when making decisions on whether to start
  // a new gogc cycle.
  if (success) {
    set_expansion_cause(cause);
    log_trace(gc)("Expanded gogc gen for %s",  gogcExpansionCause::to_string(cause));
  }
}

HeapWord* goConcurrentMarkSweepGeneration::expand_and_par_lab_allocate(gogcParGCThreadState* ps, size_t word_sz) {
  HeapWord* res = NULL;
  MutexLocker x(ParGCRareEvent_lock);
  while (true) {
    // Expansion by some other thread might make alloc OK now:
    res = ps->lab.alloc(word_sz);
    if (res != NULL) return res;
    // If there's not enough expansion space available, give up.
    if (_virtual_space.uncommitted_size() < (word_sz * HeapWordSize)) {
      return NULL;
    }
    // Otherwise, we try expansion.
    expand_for_gc_cause(word_sz*HeapWordSize, MinHeapDeltaBytes, gogcExpansionCause::_allocate_par_lab);
    // Now go around the loop and try alloc again;
    // A competing par_promote might beat us to the expansion space,
    // so we may go around the loop again if promotion fails again.
    if (GCExpandToAllocateDelayMillis > 0) {
      os::sleep(Thread::current(), GCExpandToAllocateDelayMillis, false);
    }
  }
}


bool goConcurrentMarkSweepGeneration::expand_and_ensure_spooling_space(
  PromotionInfo* promo) {
  MutexLocker x(ParGCRareEvent_lock);
  size_t refill_size_bytes = promo->refillSize() * HeapWordSize;
  while (true) {
    // Expansion by some other thread might make alloc OK now:
    if (promo->ensure_spooling_space()) {
      assert(promo->has_spooling_space(),
             "Post-condition of successful ensure_spooling_space()");
      return true;
    }
    // If there's not enough expansion space available, give up.
    if (_virtual_space.uncommitted_size() < refill_size_bytes) {
      return false;
    }
    // Otherwise, we try expansion.
    expand_for_gc_cause(refill_size_bytes, MinHeapDeltaBytes, gogcExpansionCause::_allocate_par_spooling_space);
    // Now go around the loop and try alloc again;
    // A competing allocation might beat us to the expansion space,
    // so we may go around the loop again if allocation fails again.
    if (GCExpandToAllocateDelayMillis > 0) {
      os::sleep(Thread::current(), GCExpandToAllocateDelayMillis, false);
    }
  }
}

void goConcurrentMarkSweepGeneration::shrink(size_t bytes) {
  // Only shrink if a compaction was done so that all the free space
  // in the generation is in a contiguous block at the end.
  if (did_compact()) {
    CardGeneration::shrink(bytes);
  }
}

void goConcurrentMarkSweepGeneration::assert_correct_size_change_locking() {
  assert_locked_or_safepoint(Heap_lock);
}

void goConcurrentMarkSweepGeneration::shrink_free_list_by(size_t bytes) {
  assert_locked_or_safepoint(Heap_lock);
  assert_lock_strong(freelistLock());
  log_trace(gc)("Shrinking of gogc not yet implemented");
  return;
}


// Simple ctor/dtor wrapper for accounting & timer chores around concurrent
// phases.
class gogcPhaseAccounting: public StackObj {
 public:
  gogcPhaseAccounting(gogcCollector *collector,
                     const char *title);
  ~gogcPhaseAccounting();

 private:
  gogcCollector *_collector;
  const char *_title;
  GCTraceConcTime(Info, gc) _trace_time;

 public:
  // Not MT-safe; so do not pass around these StackObj's
  // where they may be accessed by other threads.
  double wallclock_millis() {
    return TimeHelper::counter_to_millis(os::elapsed_counter() - _trace_time.start_time());
  }
};

gogcPhaseAccounting::gogcPhaseAccounting(gogcCollector *collector,
                                       const char *title) :
  _collector(collector), _title(title), _trace_time(title) {

  _collector->resetYields();
  _collector->resetTimer();
  _collector->startTimer();
  _collector->gc_timer_cm()->register_gc_concurrent_start(title);
}

gogcPhaseAccounting::~gogcPhaseAccounting() {
  _collector->gc_timer_cm()->register_gc_concurrent_end();
  _collector->stopTimer();
  log_debug(gc)("Concurrent active time: %.3fms", TimeHelper::counter_to_millis(_collector->timerTicks()));
  log_trace(gc)(" (gogc %s yielded %d times)", _title, _collector->yields());
}

// gogc work

// The common parts of gogcParInitialMarkTask and gogcParRemarkTask.
class gogcParMarkTask : public AbstractGangTask {
 protected:
  gogcCollector*     _collector;
  uint              _n_workers;
  gogcParMarkTask(const char* name, gogcCollector* collector, uint n_workers) :
      AbstractGangTask(name),
      _collector(collector),
      _n_workers(n_workers) {}
  // Work method in support of parallel rescan ... of young gen spaces
  void do_young_space_rescan(OopsInGenClosure* cl,
                             ContiguousSpace* space,
                             HeapWord** chunk_array, size_t chunk_top);
  void work_on_young_gen_roots(OopsInGenClosure* cl);
};

// Parallel initial mark task
class gogcParInitialMarkTask: public gogcParMarkTask {
  StrongRootsScope* _strong_roots_scope;
 public:
  gogcParInitialMarkTask(gogcCollector* collector, StrongRootsScope* strong_roots_scope, uint n_workers) :
      gogcParMarkTask("Scan roots and young gen for initial mark in parallel", collector, n_workers),
      _strong_roots_scope(strong_roots_scope) {}
  void work(uint worker_id);
};

// Checkpoint the roots into this generation from outside
// this generation. [Note this initial checkpoint need only
// be approximate -- we'll do a catch up phase subsequently.]
void gogcCollector::checkpointRootsInitial() {
  assert(_collectorState == InitialMarking, "Wrong collector state");
  check_correct_thread_executing();
  TracegogcMemoryManagerStats tms(_collectorState, gogcHeap::heap()->gc_cause());

  save_heap_summary();
  report_heap_summary(GCWhen::BeforeGC);

  ReferenceProcessor* rp = ref_processor();
  assert(_restart_addr == NULL, "Control point invariant");
  {
    // acquire locks for subsequent manipulations
    MutexLocker x(bitMapLock(),
                  Mutex::_no_safepoint_check_flag);
    checkpointRootsInitialWork();
    // enable ("weak") refs discovery
    rp->enable_discovery();
    _collectorState = Marking;
  }
}

void gogcCollector::checkpointRootsInitialWork() {
  assert(SafepointSynchronize::is_at_safepoint(), "world should be stopped");
  assert(_collectorState == InitialMarking, "just checking");

  // Already have locks.
  assert_lock_strong(bitMapLock());
  assert(_markBitMap.isAllClear(), "was reset at end of previous cycle");

  // Setup the verification and class unloading state for this
  // gogc collection cycle.
  setup_gogc_unloading_and_verification_state();

  GCTraceTime(Trace, gc, phases) ts("checkpointRootsInitialWork", _gc_timer_cm);

  // Reset all the PLAB chunk arrays if necessary.
  if (_survivor_plab_array != NULL && !gogcPLABRecordAlways) {
    reset_survivor_plab_arrays();
  }

  ResourceMark rm;
  HandleMark  hm;

  MarkRefsIntoClosure notOlder(_span, &_markBitMap);
  gogcHeap* heap = gogcHeap::heap();

  verify_work_stacks_empty();
  verify_overflow_empty();

  heap->ensure_parsability(false);  // fill TLABs, but no need to retire them
  // Update the saved marks which may affect the root scans.
  heap->save_marks();

  // weak reference processing has not started yet.
  ref_processor()->set_enqueuing_is_done(false);

  // Need to remember all newly created CLDs,
  // so that we can guarantee that the remark finds them.
  ClassLoaderDataGraph::remember_new_clds(true);

  // Whenever a CLD is found, it will be claimed before proceeding to mark
  // the klasses. The claimed marks need to be cleared before marking starts.
  ClassLoaderDataGraph::clear_claimed_marks();

  print_eden_and_survivor_chunk_arrays();

  {
#if COMPILER2_OR_JVMCI
    DerivedPointerTableDeactivate dpt_deact;
#endif
    if (gogcParallelInitialMarkEnabled) {
      // The parallel version.
      WorkGang* workers = heap->workers();
      assert(workers != NULL, "Need parallel worker threads.");
      uint n_workers = workers->active_workers();

      StrongRootsScope srs(n_workers);

      gogcParInitialMarkTask tsk(this, &srs, n_workers);
      initialize_sequential_subtasks_for_young_gen_rescan(n_workers);
      // If the total workers is greater than 1, then multiple workers
      // may be used at some time and the initialization has been set
      // such that the single threaded path cannot be used.
      if (workers->total_workers() > 1) {
        workers->run_task(&tsk);
      } else {
        tsk.work(0);
      }
    } else {
      // The serial version.
      CLDToOopClosure cld_closure(&notOlder, ClassLoaderData::_claim_strong);
      heap->rem_set()->prepare_for_younger_refs_iterate(false); // Not parallel.

      StrongRootsScope srs(1);

      heap->gogc_process_roots(&srs,
                             true,   // young gen as roots
                             GenCollectedHeap::ScanningOption(roots_scanning_options()),
                             should_unload_classes(),
                             &notOlder,
                             &cld_closure);
    }
  }

  // Clear mod-union table; it will be dirtied in the prologue of
  // gogc generation per each young generation collection.

  assert(_modUnionTable.isAllClear(),
       "Was cleared in most recent final checkpoint phase"
       " or no bits are set in the gc_prologue before the start of the next "
       "subsequent marking phase.");

  assert(_ct->cld_rem_set()->mod_union_is_clear(), "Must be");

  // Save the end of the used_region of the constituent generations
  // to be used to limit the extent of sweep in each generation.
  save_sweep_limits();
  verify_overflow_empty();
}

bool gogcCollector::markFromRoots() {
  // we might be tempted to assert that:
  // assert(!SafepointSynchronize::is_at_safepoint(),
  //        "inconsistent argument?");
  // However that wouldn't be right, because it's possible that
  // a safepoint is indeed in progress as a young generation
  // stop-the-world GC happens even as we mark in this generation.
  assert(_collectorState == Marking, "inconsistent state?");
  check_correct_thread_executing();
  verify_overflow_empty();

  // Weak ref discovery note: We may be discovering weak
  // refs in this generation concurrent (but interleaved) with
  // weak ref discovery by the young generation collector.

  gogcTokenSyncWithLocks ts(true, bitMapLock());
  GCTraceCPUTime tcpu;
  gogcPhaseAccounting pa(this, "Concurrent Mark");
  bool res = markFromRootsWork();
  if (res) {
    _collectorState = Precleaning;
  } else { // We failed and a foreground collection wants to take over
    assert(_foregroundGCIsActive, "internal state inconsistency");
    assert(_restart_addr == NULL,  "foreground will restart from scratch");
    log_debug(gc)("bailing out to foreground collection");
  }
  verify_overflow_empty();
  return res;
}

bool gogcCollector::markFromRootsWork() {
  // iterate over marked bits in bit map, doing a full scan and mark
  // from these roots using the following algorithm:
  // . if oop is to the right of the current scan pointer,
  //   mark corresponding bit (we'll process it later)
  // . else (oop is to left of current scan pointer)
  //   push oop on marking stack
  // . drain the marking stack

  // Note that when we do a marking step we need to hold the
  // bit map lock -- recall that direct allocation (by mutators)
  // and promotion (by the young generation collector) is also
  // marking the bit map. [the so-called allocate live policy.]
  // Because the implementation of bit map marking is not
  // robust wrt simultaneous marking of bits in the same word,
  // we need to make sure that there is no such interference
  // between concurrent such updates.

  // already have locks
  assert_lock_strong(bitMapLock());

  verify_work_stacks_empty();
  verify_overflow_empty();
  bool result = false;
  if (gogcConcurrentMTEnabled && ConcGCThreads > 0) {
    result = do_marking_mt();
  } else {
    result = do_marking_st();
  }
  return result;
}

// Forward decl
class gogcConcMarkingTask;

class gogcConcMarkingParallelTerminator: public ParallelTaskTerminator {
  gogcCollector*       _collector;
  gogcConcMarkingTask* _task;
 public:
  virtual void yield();

  // "n_threads" is the number of threads to be terminated.
  // "queue_set" is a set of work queues of other threads.
  // "collector" is the gogc collector associated with this task terminator.
  // "yield" indicates whether we need the gang as a whole to yield.
  gogcConcMarkingParallelTerminator(int n_threads, TaskQueueSetSuper* queue_set, gogcCollector* collector) :
    ParallelTaskTerminator(n_threads, queue_set),
    _collector(collector) { }

  void set_task(gogcConcMarkingTask* task) {
    _task = task;
  }
};

class gogcConcMarkingOWSTTerminator: public OWSTTaskTerminator {
  gogcCollector*       _collector;
  gogcConcMarkingTask* _task;
 public:
  virtual void yield();

  // "n_threads" is the number of threads to be terminated.
  // "queue_set" is a set of work queues of other threads.
  // "collector" is the gogc collector associated with this task terminator.
  // "yield" indicates whether we need the gang as a whole to yield.
  gogcConcMarkingOWSTTerminator(int n_threads, TaskQueueSetSuper* queue_set, gogcCollector* collector) :
    OWSTTaskTerminator(n_threads, queue_set),
    _collector(collector) { }

  void set_task(gogcConcMarkingTask* task) {
    _task = task;
  }
};

class gogcConcMarkingTaskTerminator {
 private:
  ParallelTaskTerminator* _term;
 public:
  gogcConcMarkingTaskTerminator(int n_threads, TaskQueueSetSuper* queue_set, gogcCollector* collector) {
    if (UseOWSTTaskTerminator) {
      _term = new gogcConcMarkingOWSTTerminator(n_threads, queue_set, collector);
    } else {
      _term = new gogcConcMarkingParallelTerminator(n_threads, queue_set, collector);
    }
  }
  ~gogcConcMarkingTaskTerminator() {
    assert(_term != NULL, "Must not be NULL");
    delete _term;
  }

  void set_task(gogcConcMarkingTask* task);
  ParallelTaskTerminator* terminator() const { return _term; }
};

class gogcConcMarkingTerminatorTerminator: public TerminatorTerminator {
  gogcConcMarkingTask* _task;
 public:
  bool should_exit_termination();
  void set_task(gogcConcMarkingTask* task) {
    _task = task;
  }
};

// MT Concurrent Marking Task
class gogcConcMarkingTask: public YieldingFlexibleGangTask {
  gogcCollector*             _collector;
  uint                      _n_workers;      // requested/desired # workers
  bool                      _result;
  CompactibleFreeListSpace* _gogc_space;
  char                      _pad_front[64];   // padding to ...
  HeapWord* volatile        _global_finger;   // ... avoid sharing cache line
  char                      _pad_back[64];
  HeapWord*                 _restart_addr;

  //  Exposed here for yielding support
  Mutex* const _bit_map_lock;

  // The per thread work queues, available here for stealing
  OopTaskQueueSet*  _task_queues;

  // Termination (and yielding) support
  gogcConcMarkingTaskTerminator       _term;
  gogcConcMarkingTerminatorTerminator _term_term;

 public:
  gogcConcMarkingTask(gogcCollector* collector,
                 CompactibleFreeListSpace* gogc_space,
                 YieldingFlexibleWorkGang* workers,
                 OopTaskQueueSet* task_queues):
    YieldingFlexibleGangTask("Concurrent marking done multi-threaded"),
    _collector(collector),
    _n_workers(0),
    _result(true),
    _gogc_space(gogc_space),
    _bit_map_lock(collector->bitMapLock()),
    _task_queues(task_queues),
    _term(_n_workers, task_queues, _collector)
  {
    _requested_size = _n_workers;
    _term.set_task(this);
    _term_term.set_task(this);
    _restart_addr = _global_finger = _gogc_space->bottom();
  }


  OopTaskQueueSet* task_queues()  { return _task_queues; }

  OopTaskQueue* work_queue(int i) { return task_queues()->queue(i); }

  HeapWord* volatile* global_finger_addr() { return &_global_finger; }

  ParallelTaskTerminator* terminator() { return _term.terminator(); }

  virtual void set_for_termination(uint active_workers) {
    terminator()->reset_for_reuse(active_workers);
  }

  void work(uint worker_id);
  bool should_yield() {
    return    ConcurrentMarkSweepThread::should_yield()
           && !_collector->foregroundGCIsActive();
  }

  virtual void coordinator_yield();  // stuff done by coordinator
  bool result() { return _result; }

  void reset(HeapWord* ra) {
    assert(_global_finger >= _gogc_space->end(),  "Postcondition of ::work(i)");
    _restart_addr = _global_finger = ra;
    _term.terminator()->reset_for_reuse();
  }

  static bool get_work_from_overflow_stack(gogcMarkStack* ovflw_stk,
                                           OopTaskQueue* work_q);

 private:
  void do_scan_and_mark(int i, CompactibleFreeListSpace* sp);
  void do_work_steal(int i);
  void bump_global_finger(HeapWord* f);
};

bool gogcConcMarkingTerminatorTerminator::should_exit_termination() {
  assert(_task != NULL, "Error");
  return _task->yielding();
  // Note that we do not need the disjunct || _task->should_yield() above
  // because we want terminating threads to yield only if the task
  // is already in the midst of yielding, which happens only after at least one
  // thread has yielded.
}

void gogcConcMarkingParallelTerminator::yield() {
  if (_task->should_yield()) {
    _task->yield();
  } else {
    ParallelTaskTerminator::yield();
  }
}

void gogcConcMarkingOWSTTerminator::yield() {
  if (_task->should_yield()) {
    _task->yield();
  } else {
    OWSTTaskTerminator::yield();
  }
}

void gogcConcMarkingTaskTerminator::set_task(gogcConcMarkingTask* task) {
  if (UseOWSTTaskTerminator) {
    ((gogcConcMarkingOWSTTerminator*)_term)->set_task(task);
  } else {
    ((gogcConcMarkingParallelTerminator*)_term)->set_task(task);
  }
}

////////////////////////////////////////////////////////////////
// Concurrent Marking Algorithm Sketch
////////////////////////////////////////////////////////////////
// Until all tasks exhausted (both spaces):
// -- claim next available chunk
// -- bump global finger via CAS
// -- find first object that starts in this chunk
//    and start scanning bitmap from that position
// -- scan marked objects for oops
// -- CAS-mark target, and if successful:
//    . if target oop is above global finger (volatile read)
//      nothing to do
//    . if target oop is in chunk and above local finger
//        then nothing to do
//    . else push on work-queue
// -- Deal with possible overflow issues:
//    . local work-queue overflow causes stuff to be pushed on
//      global (common) overflow queue
//    . always first empty local work queue
//    . then get a batch of oops from global work queue if any
//    . then do work stealing
// -- When all tasks claimed (both spaces)
//    and local work queue empty,
//    then in a loop do:
//    . check global overflow stack; steal a batch of oops and trace
//    . try to steal from other threads oif GOS is empty
//    . if neither is available, offer termination
// -- Terminate and return result
//
void gogcConcMarkingTask::work(uint worker_id) {
  elapsedTimer _timer;
  ResourceMark rm;
  HandleMark hm;

  DEBUG_ONLY(_collector->verify_overflow_empty();)

  // Before we begin work, our work queue should be empty
  assert(work_queue(worker_id)->size() == 0, "Expected to be empty");
  // Scan the bitmap covering _gogc_space, tracing through grey objects.
  _timer.start();
  do_scan_and_mark(worker_id, _gogc_space);
  _timer.stop();
  log_trace(gc, task)("Finished gogc space scanning in %dth thread: %3.3f sec", worker_id, _timer.seconds());

  // ... do work stealing
  _timer.reset();
  _timer.start();
  do_work_steal(worker_id);
  _timer.stop();
  log_trace(gc, task)("Finished work stealing in %dth thread: %3.3f sec", worker_id, _timer.seconds());
  assert(_collector->_markStack.isEmpty(), "Should have been emptied");
  assert(work_queue(worker_id)->size() == 0, "Should have been emptied");
  // Note that under the current task protocol, the
  // following assertion is true even of the spaces
  // expanded since the completion of the concurrent
  // marking. XXX This will likely change under a strict
  // ABORT semantics.
  // After perm removal the comparison was changed to
  // greater than or equal to from strictly greater than.
  // Before perm removal the highest address sweep would
  // have been at the end of perm gen but now is at the
  // end of the tenured gen.
  assert(_global_finger >=  _gogc_space->end(),
         "All tasks have been completed");
  DEBUG_ONLY(_collector->verify_overflow_empty();)
}

void gogcConcMarkingTask::bump_global_finger(HeapWord* f) {
  HeapWord* read = _global_finger;
  HeapWord* cur  = read;
  while (f > read) {
    cur = read;
    read = Atomic::cmpxchg(f, &_global_finger, cur);
    if (cur == read) {
      // our cas succeeded
      assert(_global_finger >= f, "protocol consistency");
      break;
    }
  }
}

// This is really inefficient, and should be redone by
// using (not yet available) block-read and -write interfaces to the
// stack and the work_queue. XXX FIX ME !!!
bool gogcConcMarkingTask::get_work_from_overflow_stack(gogcMarkStack* ovflw_stk,
                                                      OopTaskQueue* work_q) {
  // Fast lock-free check
  if (ovflw_stk->length() == 0) {
    return false;
  }
  assert(work_q->size() == 0, "Shouldn't steal");
  MutexLocker ml(ovflw_stk->par_lock(),
                 Mutex::_no_safepoint_check_flag);
  // Grab up to 1/4 the size of the work queue
  size_t num = MIN2((size_t)(work_q->max_elems() - work_q->size())/4,
                    (size_t)ParGCDesiredObjsFromOverflowList);
  num = MIN2(num, ovflw_stk->length());
  for (int i = (int) num; i > 0; i--) {
    oop cur = ovflw_stk->pop();
    assert(cur != NULL, "Counted wrong?");
    work_q->push(cur);
  }
  return num > 0;
}

void gogcConcMarkingTask::do_scan_and_mark(int i, CompactibleFreeListSpace* sp) {
  SequentialSubTasksDone* pst = sp->conc_par_seq_tasks();
  int n_tasks = pst->n_tasks();
  // We allow that there may be no tasks to do here because
  // we are restarting after a stack overflow.
  assert(pst->valid() || n_tasks == 0, "Uninitialized use?");
  uint nth_task = 0;

  HeapWord* aligned_start = sp->bottom();
  if (sp->used_region().contains(_restart_addr)) {
    // Align down to a card boundary for the start of 0th task
    // for this space.
    aligned_start = align_down(_restart_addr, CardTable::card_size);
  }

  size_t chunk_size = sp->marking_task_size();
  while (pst->try_claim_task(/* reference */ nth_task)) {
    // Having claimed the nth task in this space,
    // compute the chunk that it corresponds to:
    MemRegion span = MemRegion(aligned_start + nth_task*chunk_size,
                               aligned_start + (nth_task+1)*chunk_size);
    // Try and bump the global finger via a CAS;
    // note that we need to do the global finger bump
    // _before_ taking the intersection below, because
    // the task corresponding to that region will be
    // deemed done even if the used_region() expands
    // because of allocation -- as it almost certainly will
    // during start-up while the threads yield in the
    // closure below.
    HeapWord* finger = span.end();
    bump_global_finger(finger);   // atomically
    // There are null tasks here corresponding to chunks
    // beyond the "top" address of the space.
    span = span.intersection(sp->used_region());
    if (!span.is_empty()) {  // Non-null task
      HeapWord* prev_obj;
      assert(!span.contains(_restart_addr) || nth_task == 0,
             "Inconsistency");
      if (nth_task == 0) {
        // For the 0th task, we'll not need to compute a block_start.
        if (span.contains(_restart_addr)) {
          // In the case of a restart because of stack overflow,
          // we might additionally skip a chunk prefix.
          prev_obj = _restart_addr;
        } else {
          prev_obj = span.start();
        }
      } else {
        // We want to skip the first object because
        // the protocol is to scan any object in its entirety
        // that _starts_ in this span; a fortiori, any
        // object starting in an earlier span is scanned
        // as part of an earlier claimed task.
        // Below we use the "careful" version of block_start
        // so we do not try to navigate uninitialized objects.
        prev_obj = sp->block_start_careful(span.start());
        // Below we use a variant of block_size that uses the
        // Printezis bits to avoid waiting for allocated
        // objects to become initialized/parsable.
        while (prev_obj < span.start()) {
          size_t sz = sp->block_size_no_stall(prev_obj, _collector);
          if (sz > 0) {
            prev_obj += sz;
          } else {
            // In this case we may end up doing a bit of redundant
            // scanning, but that appears unavoidable, short of
            // locking the free list locks; see bug 6324141.
            break;
          }
        }
      }
      if (prev_obj < span.end()) {
        MemRegion my_span = MemRegion(prev_obj, span.end());
        // Do the marking work within a non-empty span --
        // the last argument to the constructor indicates whether the
        // iteration should be incremental with periodic yields.
        ParMarkFromRootsClosure cl(this, _collector, my_span,
                                   &_collector->_markBitMap,
                                   work_queue(i),
                                   &_collector->_markStack);
        _collector->_markBitMap.iterate(&cl, my_span.start(), my_span.end());
      } // else nothing to do for this task
    }   // else nothing to do for this task
  }
  // We'd be tempted to assert here that since there are no
  // more tasks left to claim in this space, the global_finger
  // must exceed space->top() and a fortiori space->end(). However,
  // that would not quite be correct because the bumping of
  // global_finger occurs strictly after the claiming of a task,
  // so by the time we reach here the global finger may not yet
  // have been bumped up by the thread that claimed the last
  // task.
  pst->all_tasks_completed();
}

class ParConcMarkingClosure: public MetadataVisitingOopIterateClosure {
 private:
  gogcCollector* _collector;
  gogcConcMarkingTask* _task;
  MemRegion     _span;
  gogcBitMap*    _bit_map;
  gogcMarkStack* _overflow_stack;
  OopTaskQueue* _work_queue;
 protected:
  DO_OOP_WORK_DEFN
 public:
  ParConcMarkingClosure(gogcCollector* collector, gogcConcMarkingTask* task, OopTaskQueue* work_queue,
                        gogcBitMap* bit_map, gogcMarkStack* overflow_stack):
    MetadataVisitingOopIterateClosure(collector->ref_processor()),
    _collector(collector),
    _task(task),
    _span(collector->_span),
    _bit_map(bit_map),
    _overflow_stack(overflow_stack),
    _work_queue(work_queue)
  { }
  virtual void do_oop(oop* p);
  virtual void do_oop(narrowOop* p);

  void trim_queue(size_t max);
  void handle_stack_overflow(HeapWord* lost);
  void do_yield_check() {
    if (_task->should_yield()) {
      _task->yield();
    }
  }
};

DO_OOP_WORK_IMPL(ParConcMarkingClosure)

// Grey object scanning during work stealing phase --
// the salient assumption here is that any references
// that are in these stolen objects being scanned must
// already have been initialized (else they would not have
// been published), so we do not need to check for
// uninitialized objects before pushing here.
void ParConcMarkingClosure::do_oop(oop obj) {
  assert(oopDesc::is_oop_or_null(obj, true), "Expected an oop or NULL at " PTR_FORMAT, p2i(obj));
  HeapWord* addr = (HeapWord*)obj;
  // Check if oop points into the gogc generation
  // and is not marked
  if (_span.contains(addr) && !_bit_map->isMarked(addr)) {
    // a white object ...
    // If we manage to "claim" the object, by being the
    // first thread to mark it, then we push it on our
    // marking stack
    if (_bit_map->par_mark(addr)) {     // ... now grey
      // push on work queue (grey set)
      bool simulate_overflow = false;
      NOT_PRODUCT(
        if (gogcMarkStackOverflowALot &&
            _collector->simulate_overflow()) {
          // simulate a stack overflow
          simulate_overflow = true;
        }
      )
      if (simulate_overflow ||
          !(_work_queue->push(obj) || _overflow_stack->par_push(obj))) {
        // stack overflow
        log_trace(gc)("gogc marking stack overflow (benign) at " SIZE_FORMAT, _overflow_stack->capacity());
        // We cannot assert that the overflow stack is full because
        // it may have been emptied since.
        assert(simulate_overflow ||
               _work_queue->size() == _work_queue->max_elems(),
              "Else push should have succeeded");
        handle_stack_overflow(addr);
      }
    } // Else, some other thread got there first
    do_yield_check();
  }
}
