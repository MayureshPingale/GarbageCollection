#ifndef SHARE_GC_gogc_CONCURRENTMARKSWEEPGENERATION_INLINE_HPP
#define SHARE_GC_gogc_CONCURRENTMARKSWEEPGENERATION_INLINE_HPP

#include "gc/gogc/gogcHeap.hpp"
#include "gc/gogc/gogcLockVerifier.hpp"
#include "gc/gogc/compactibleFreeListSpace.inline.hpp"
#include "gc/gogc/concurrentMarkSweepGeneration.hpp"
#include "gc/gogc/goConcurrentMarkSweepThread.hpp"
#include "gc/gogc/parNewGeneration.hpp"
#include "gc/shared/gcUtil.hpp"
#include "utilities/align.hpp"
#include "utilities/bitMap.inline.hpp"

inline void gogcBitMap::clear_all() {
  assert_locked();
  _bm.clear_large();
  return;
}

inline size_t gogcBitMap::heapWordToOffset(HeapWord* addr) const {
  return (pointer_delta(addr, _bmStartWord)) >> _shifter;
}

inline HeapWord* gogcBitMap::offsetToHeapWord(size_t offset) const {
  return _bmStartWord + (offset << _shifter);
}

inline size_t gogcBitMap::heapWordDiffToOffsetDiff(size_t diff) const {
  assert((diff & ((1 << _shifter) - 1)) == 0, "argument check");
  return diff >> _shifter;
}

inline void gogcBitMap::mark(HeapWord* addr) {
  assert_locked();
  assert(_bmStartWord <= addr && addr < (_bmStartWord + _bmWordSize),
         "outside underlying space?");
  _bm.set_bit(heapWordToOffset(addr));
}

inline bool gogcBitMap::par_mark(HeapWord* addr) {
  assert_locked();
  assert(_bmStartWord <= addr && addr < (_bmStartWord + _bmWordSize),
         "outside underlying space?");
  return _bm.par_at_put(heapWordToOffset(addr), true);
}

inline void gogcBitMap::par_clear(HeapWord* addr) {
  assert_locked();
  assert(_bmStartWord <= addr && addr < (_bmStartWord + _bmWordSize),
         "outside underlying space?");
  _bm.par_at_put(heapWordToOffset(addr), false);
}

inline void gogcBitMap::mark_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  // Range size is usually just 1 bit.
  _bm.set_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                BitMap::small_range);
}

inline void gogcBitMap::clear_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  // Range size is usually just 1 bit.
  _bm.clear_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                  BitMap::small_range);
}

inline void gogcBitMap::par_mark_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  _bm.par_set_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                    BitMap::small_range);
}

inline void gogcBitMap::par_clear_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  _bm.par_clear_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                      BitMap::small_range);
}

inline void gogcBitMap::mark_large_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  _bm.set_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                BitMap::large_range);
}

inline void gogcBitMap::clear_large_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  _bm.clear_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                  BitMap::large_range);
}

inline void gogcBitMap::par_mark_large_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  _bm.par_set_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                    BitMap::large_range);
}

inline void gogcBitMap::par_clear_large_range(MemRegion mr) {
  NOT_PRODUCT(region_invariant(mr));
  _bm.par_clear_range(heapWordToOffset(mr.start()), heapWordToOffset(mr.end()),
                      BitMap::large_range);
}


inline MemRegion gogcBitMap::getAndClearMarkedRegion(HeapWord* addr) {
  return getAndClearMarkedRegion(addr, endWord());
}


inline MemRegion gogcBitMap::getAndClearMarkedRegion(HeapWord* start_addr,
                                                    HeapWord* end_addr) {
  HeapWord *start, *end;
  assert_locked();
  start = getNextMarkedWordAddress  (start_addr, end_addr);
  end   = getNextUnmarkedWordAddress(start,      end_addr);
  assert(start <= end, "Consistency check");
  MemRegion mr(start, end);
  if (!mr.is_empty()) {
    clear_range(mr);
  }
  return mr;
}

inline bool gogcBitMap::isMarked(HeapWord* addr) const {
  assert_locked();
  assert(_bmStartWord <= addr && addr < (_bmStartWord + _bmWordSize),
         "outside underlying space?");
  return _bm.at(heapWordToOffset(addr));
}


inline bool gogcBitMap::par_isMarked(HeapWord* addr) const {
  assert(_bmStartWord <= addr && addr < (_bmStartWord + _bmWordSize),
         "outside underlying space?");
  return _bm.at(heapWordToOffset(addr));
}


inline bool gogcBitMap::isUnmarked(HeapWord* addr) const {
  assert_locked();
  assert(_bmStartWord <= addr && addr < (_bmStartWord + _bmWordSize),
         "outside underlying space?");
  return !_bm.at(heapWordToOffset(addr));
}

inline HeapWord* gogcBitMap::getNextMarkedWordAddress(HeapWord* addr) const {
  return getNextMarkedWordAddress(addr, endWord());
}

inline HeapWord* gogcBitMap::getNextMarkedWordAddress(
  HeapWord* start_addr, HeapWord* end_addr) const {
  assert_locked();
  size_t nextOffset = _bm.get_next_one_offset(
                        heapWordToOffset(start_addr),
                        heapWordToOffset(end_addr));
  HeapWord* nextAddr = offsetToHeapWord(nextOffset);
  assert(nextAddr >= start_addr &&
         nextAddr <= end_addr, "get_next_one postcondition");
  assert((nextAddr == end_addr) ||
         isMarked(nextAddr), "get_next_one postcondition");
  return nextAddr;
}


inline HeapWord* gogcBitMap::getNextUnmarkedWordAddress(HeapWord* addr) const {
  return getNextUnmarkedWordAddress(addr, endWord());
}


inline HeapWord* gogcBitMap::getNextUnmarkedWordAddress(
  HeapWord* start_addr, HeapWord* end_addr) const {
  assert_locked();
  size_t nextOffset = _bm.get_next_zero_offset(
                        heapWordToOffset(start_addr),
                        heapWordToOffset(end_addr));
  HeapWord* nextAddr = offsetToHeapWord(nextOffset);
  assert(nextAddr >= start_addr &&
         nextAddr <= end_addr, "get_next_zero postcondition");
  assert((nextAddr == end_addr) ||
          isUnmarked(nextAddr), "get_next_zero postcondition");
  return nextAddr;
}

inline bool gogcBitMap::isAllClear() const {
  assert_locked();
  return getNextMarkedWordAddress(startWord()) >= endWord();
}

inline void gogcBitMap::iterate(BitMapClosure* cl, HeapWord* left,
                            HeapWord* right) {
  assert_locked();
  left = MAX2(_bmStartWord, left);
  right = MIN2(_bmStartWord + _bmWordSize, right);
  if (right > left) {
    _bm.iterate(cl, heapWordToOffset(left), heapWordToOffset(right));
  }
}

inline void gogcCollector::save_sweep_limits() {
  _gogcGen->save_sweep_limit();
}

inline bool gogcCollector::is_dead_obj(oop obj) const {
  HeapWord* addr = (HeapWord*)obj;
  assert((_gogcGen->gogcSpace()->is_in_reserved(addr)
          && _gogcGen->gogcSpace()->block_is_obj(addr)),
         "must be object");
  return  should_unload_classes() &&
          _collectorState == Sweeping &&
         !_markBitMap.isMarked(addr);
}

inline bool gogcCollector::should_abort_preclean() const {
  return _collectorState == AbortablePreclean &&
         (_abort_preclean || _foregroundGCIsActive ||
          gogcHeap::heap()->incremental_collection_will_fail(true /* consult_young */));
}

inline size_t gogcCollector::get_eden_used() const {
  return _young_gen->eden()->used();
}

inline size_t gogcCollector::get_eden_capacity() const {
  return _young_gen->eden()->capacity();
}

inline bool gogcStats::valid() const {
  return _valid_bits == _ALL_VALID;
}

inline void gogcStats::record_gc0_begin() {
  if (_gc0_begin_time.is_updated()) {
    float last_gc0_period = _gc0_begin_time.seconds();
    _gc0_period = AdaptiveWeightedAverage::exp_avg(_gc0_period,
      last_gc0_period, _gc0_alpha);
    _gc0_alpha = _saved_alpha;
    _valid_bits |= _GC0_VALID;
  }
  _gogc_used_at_gc0_begin = _gogc_gen->gogcSpace()->used();

  _gc0_begin_time.update();
}

inline void gogcStats::record_gc0_end(size_t gogc_gen_bytes_used) {
  float last_gc0_duration = _gc0_begin_time.seconds();
  _gc0_duration = AdaptiveWeightedAverage::exp_avg(_gc0_duration,
    last_gc0_duration, _gc0_alpha);

  // Amount promoted.
  _gogc_used_at_gc0_end = gogc_gen_bytes_used;

  size_t promoted_bytes = 0;
  if (_gogc_used_at_gc0_end >= _gogc_used_at_gc0_begin) {
    promoted_bytes = _gogc_used_at_gc0_end - _gogc_used_at_gc0_begin;
  }

  _gogc_gen->gc_stats()->avg_promoted()->sample(promoted_bytes);
  _gc0_promoted = (size_t) _gogc_gen->gc_stats()->avg_promoted()->average();

  // Amount directly allocated.
  size_t allocated_bytes = _gogc_gen->direct_allocated_words() * HeapWordSize;
  _gogc_gen->reset_direct_allocated_words();
  _gogc_allocated = AdaptiveWeightedAverage::exp_avg(_gogc_allocated,
    allocated_bytes, _gc0_alpha);
}

inline void gogcStats::record_gogc_begin() {
  _gogc_timer.stop();

  _gogc_used_at_gogc_begin = _gogc_used_at_gc0_end;

  _gogc_period = AdaptiveWeightedAverage::exp_avg((float)_gogc_period,
    (float) _gogc_timer.seconds(), _gogc_alpha);
  _gogc_begin_time.update();

  _gogc_timer.reset();
  _gogc_timer.start();
}

inline void gogcStats::record_gogc_end() {
  _gogc_timer.stop();

  float cur_duration = _gogc_timer.seconds();
  _gogc_duration = AdaptiveWeightedAverage::exp_avg(_gogc_duration,
    cur_duration, _gogc_alpha);

  _gogc_end_time.update();
  _gogc_alpha = _saved_alpha;
  _allow_duty_cycle_reduction = true;
  _valid_bits |= _gogc_VALID;

  _gogc_timer.start();
}

inline double gogcStats::gogc_time_since_begin() const {
  return _gogc_begin_time.seconds();
}

inline double gogcStats::gogc_time_since_end() const {
  return _gogc_end_time.seconds();
}

inline double gogcStats::promotion_rate() const {
  assert(valid(), "statistics not valid yet");
  return gc0_promoted() / gc0_period();
}

inline double gogcStats::gogc_allocation_rate() const {
  assert(valid(), "statistics not valid yet");
  return gogc_allocated() / gc0_period();
}

inline double gogcStats::gogc_consumption_rate() const {
  assert(valid(), "statistics not valid yet");
  return (gc0_promoted() + gogc_allocated()) / gc0_period();
}

inline void ConcurrentMarkSweepGeneration::save_sweep_limit() {
  gogcSpace()->save_sweep_limit();
}

inline MemRegion ConcurrentMarkSweepGeneration::used_region_at_save_marks() const {
  return _gogcSpace->used_region_at_save_marks();
}

template <typename OopClosureType>
void ConcurrentMarkSweepGeneration::oop_since_save_marks_iterate(OopClosureType* cl) {
  cl->set_generation(this);
  gogcSpace()->oop_since_save_marks_iterate(cl);
  cl->reset_generation();
  save_marks();
}

inline void MarkFromRootsClosure::do_yield_check() {
  if (goConcurrentMarkSweepThread::should_yield() &&
      !_collector->foregroundGCIsActive() &&
      _yield) {
    do_yield_work();
  }
}

inline void ParMarkFromRootsClosure::do_yield_check() {
  if (goConcurrentMarkSweepThread::should_yield() &&
      !_collector->foregroundGCIsActive()) {
    do_yield_work();
  }
}

inline void PushOrMarkClosure::do_yield_check() {
  _parent->do_yield_check();
}

inline void ParPushOrMarkClosure::do_yield_check() {
  _parent->do_yield_check();
}


inline bool ScanMarkedObjectsAgainCarefullyClosure::do_yield_check() {
  if (goConcurrentMarkSweepThread::should_yield() &&
      !_collector->foregroundGCIsActive() &&
      _yield) {

    _collector->sample_eden();
    do_yield_work();
    _collector->sample_eden();
    return _collector->should_abort_preclean();
  }
  return false;
}

inline void SurvivorSpacePrecleanClosure::do_yield_check() {
  if (goConcurrentMarkSweepThread::should_yield() &&
      !_collector->foregroundGCIsActive() &&
      _yield) {
    _collector->sample_eden();
    do_yield_work();
    _collector->sample_eden();
  }
}

inline void SweepClosure::do_yield_check(HeapWord* addr) {
  if (goConcurrentMarkSweepThread::should_yield() &&
      !_collector->foregroundGCIsActive() &&
      _yield) {
    do_yield_work(addr);
  }
}

inline void MarkRefsIntoAndScanClosure::do_yield_check() {
  
  if (_yield &&
      !_collector->foregroundGCIsActive() &&
      goConcurrentMarkSweepThread::should_yield()) {
    do_yield_work();
  }
}


inline void ModUnionClosure::do_MemRegion(MemRegion mr) {
  MemRegion mr2(mr.start(), align_up(mr.end(),
                CardTable::card_size /* bytes */));
  _t->mark_range(mr2);
}

inline void ModUnionClosurePar::do_MemRegion(MemRegion mr) {
  MemRegion mr2(mr.start(), align_up(mr.end(),
                CardTable::card_size /* bytes */));
  _t->par_mark_range(mr2);
}

#endif // SHARE_GC_gogc_CONCURRENTMARKSWEEPGENERATION_INLINE_HPP
