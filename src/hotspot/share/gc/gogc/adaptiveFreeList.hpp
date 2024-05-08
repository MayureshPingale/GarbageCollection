#ifndef SHARE_GC_CMS_ADAPTIVEFREELIST_HPP
#define SHARE_GC_CMS_ADAPTIVEFREELIST_HPP

#include "gc/cms/allocationStats.hpp"
#include "memory/freeList.hpp"

class CompactibleFreeListSpace;
class Mutex;

template <class Chunk>
class AdaptiveFreeList : public FreeList<Chunk> {
  friend class CompactibleFreeListSpace;
  friend class VMStructs;

  size_t        _hint;          
  AllocationStats _allocation_stats; 
 public:

  AdaptiveFreeList();

  using FreeList<Chunk>::assert_proper_lock_protection;
#ifdef ASSERT
  using FreeList<Chunk>::protecting_lock;
#endif
  using FreeList<Chunk>::count;
  using FreeList<Chunk>::size;
  using FreeList<Chunk>::verify_chunk_in_free_list;
  using FreeList<Chunk>::getFirstNChunksFromList;
  using FreeList<Chunk>::print_on;
  void return_chunk_at_head(Chunk* fc, bool record_return);
  void return_chunk_at_head(Chunk* fc);
  void return_chunk_at_tail(Chunk* fc, bool record_return);
  void return_chunk_at_tail(Chunk* fc);
  using FreeList<Chunk>::return_chunk_at_tail;
  using FreeList<Chunk>::remove_chunk;
  using FreeList<Chunk>::prepend;
  using FreeList<Chunk>::print_labels_on;
  using FreeList<Chunk>::get_chunk_at_head;

  void initialize();
 void reset(size_t hint);

  void print_on(outputStream* st, const char* c = NULL) const;

  size_t hint() const {
    return _hint;
  }
  void set_hint(size_t v) {
    assert_proper_lock_protection();
    assert(v == 0 || size() < v, "Bad hint");
    _hint = v;
  }

  size_t get_better_size();
 void init_statistics(bool split_birth = false);

  AllocationStats* allocation_stats() {
    assert_proper_lock_protection();
    return &_allocation_stats;
  }

  ssize_t desired() const {
    return _allocation_stats.desired();
  }
  void set_desired(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_desired(v);
  }
  void compute_desired(float inter_sweep_current,
                       float inter_sweep_estimate,
                       float intra_sweep_estimate) {
    assert_proper_lock_protection();
    _allocation_stats.compute_desired(count(),
                                      inter_sweep_current,
                                      inter_sweep_estimate,
                                      intra_sweep_estimate);
  }
  ssize_t coal_desired() const {
    return _allocation_stats.coal_desired();
  }
  void set_coal_desired(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_coal_desired(v);
  }

  ssize_t surplus() const {
    return _allocation_stats.surplus();
  }
  void set_surplus(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_surplus(v);
  }
  void increment_surplus() {
    assert_proper_lock_protection();
    _allocation_stats.increment_surplus();
  }
  void decrement_surplus() {
    assert_proper_lock_protection();
    _allocation_stats.decrement_surplus();
  }

  ssize_t bfr_surp() const {
    return _allocation_stats.bfr_surp();
  }
  void set_bfr_surp(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_bfr_surp(v);
  }
  ssize_t prev_sweep() const {
    return _allocation_stats.prev_sweep();
  }
  void set_prev_sweep(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_prev_sweep(v);
  }
  ssize_t before_sweep() const {
    return _allocation_stats.before_sweep();
  }
  void set_before_sweep(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_before_sweep(v);
  }

  ssize_t coal_births() const {
    return _allocation_stats.coal_births();
  }
  void set_coal_births(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_coal_births(v);
  }
  void increment_coal_births() {
    assert_proper_lock_protection();
    _allocation_stats.increment_coal_births();
  }

  ssize_t coal_deaths() const {
    return _allocation_stats.coal_deaths();
  }
  void set_coal_deaths(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_coal_deaths(v);
  }
  void increment_coal_deaths() {
    assert_proper_lock_protection();
    _allocation_stats.increment_coal_deaths();
  }

  ssize_t split_births() const {
    return _allocation_stats.split_births();
  }
  void set_split_births(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_split_births(v);
  }
  void increment_split_births() {
    assert_proper_lock_protection();
    _allocation_stats.increment_split_births();
  }

  ssize_t split_deaths() const {
    return _allocation_stats.split_deaths();
  }
  void set_split_deaths(ssize_t v) {
    assert_proper_lock_protection();
    _allocation_stats.set_split_deaths(v);
  }
  void increment_split_deaths() {
    assert_proper_lock_protection();
    _allocation_stats.increment_split_deaths();
  }

#ifndef PRODUCT
  size_t returned_bytes() const { return _allocation_stats.returned_bytes(); }
  void set_returned_bytes(size_t v) { _allocation_stats.set_returned_bytes(v); }
  void increment_returned_bytes_by(size_t v) {
    _allocation_stats.set_returned_bytes(_allocation_stats.returned_bytes() + v);
  }
  void verify_stats() const;
#endif  
};

#endif 
