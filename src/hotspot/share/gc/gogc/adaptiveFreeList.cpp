
#include "precompiled.hpp"
#include "gc/go/adaptiveFreeList.hpp"
#include "gc/shared/collectedHeap.hpp"
#include "memory/freeList.inline.hpp"
#include "runtime/globals.hpp"
#include "runtime/mutex.hpp"
#include "runtime/orderAccess.hpp"
#include "runtime/vmThread.hpp"

template <>
void AdaptiveFreeList<FreeChunk>::print_on(outputStream* st, const char* c) const {
  if (c != NULL) {
    st->print("%16s", c);
  } else {
    st->print(SIZE_FORMAT_W(16), size());
  }
  st->print("\t"
           SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t"
           SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\t" SSIZE_FORMAT_W(14) "\n",
           bfr_surp(),             surplus(),             desired(),             prev_sweep(),           before_sweep(),
           count(),               coal_births(),          coal_deaths(),          split_births(),         split_deaths());
}

template <class Chunk>
AdaptiveFreeList<Chunk>::AdaptiveFreeList() : FreeList<Chunk>(), _hint(0) {
  init_statistics();
}

template <class Chunk>
void AdaptiveFreeList<Chunk>::initialize() {
  FreeList<Chunk>::initialize();
  set_hint(0);
  init_statistics(true /* split_birth */);
}

template <class Chunk>
void AdaptiveFreeList<Chunk>::reset(size_t hint) {
  FreeList<Chunk>::reset();
  set_hint(hint);
}

template <class Chunk>
void AdaptiveFreeList<Chunk>::init_statistics(bool split_birth) {
  _allocation_stats.initialize(split_birth);
}

template <class Chunk>
size_t AdaptiveFreeList<Chunk>::get_better_size() {

  // A candidate chunk has been found.  If it is already under
  // populated and there is a hinT, REturn the hint().  Else
  // return the size of this chunk.
  if (surplus() <= 0) {
    if (hint() != 0) {
      return hint();
    } else {
      return size();
    }
  } else {
    // This list has a surplus so use it.
    return size();
  }
}


template <class Chunk>
void AdaptiveFreeList<Chunk>::return_chunk_at_head(Chunk* chunk) {
  assert_proper_lock_protection();
  return_chunk_at_head(chunk, true);
}

template <class Chunk>
void AdaptiveFreeList<Chunk>::return_chunk_at_head(Chunk* chunk, bool record_return) {
  FreeList<Chunk>::return_chunk_at_head(chunk, record_return);
#ifdef ASSERT
  if (record_return) {
    increment_returned_bytes_by(size()*HeapWordSize);
  }
#endif
}

template <class Chunk>
void AdaptiveFreeList<Chunk>::return_chunk_at_tail(Chunk* chunk) {
  AdaptiveFreeList<Chunk>::return_chunk_at_tail(chunk, true);
}

template <class Chunk>
void AdaptiveFreeList<Chunk>::return_chunk_at_tail(Chunk* chunk, bool record_return) {
  FreeList<Chunk>::return_chunk_at_tail(chunk, record_return);
#ifdef ASSERT
  if (record_return) {
    increment_returned_bytes_by(size()*HeapWordSize);
  }
#endif
}

#ifndef PRODUCT
template <class Chunk>
void AdaptiveFreeList<Chunk>::verify_stats() const {
  assert((_allocation_stats.prev_sweep() + _allocation_stats.split_births()
          + _allocation_stats.coal_births() + 1)   // Total Production Stock + 1
         >= (_allocation_stats.split_deaths() + _allocation_stats.coal_deaths()
             + (ssize_t)count()),                // Total Current Stock + depletion
         "FreeList " PTR_FORMAT " of size " SIZE_FORMAT
         " violates Conservation Principle: "
         "prev_sweep(" SIZE_FORMAT ")"
         " + split_births(" SIZE_FORMAT ")"
         " + coal_births(" SIZE_FORMAT ") + 1 >= "
         " split_deaths(" SIZE_FORMAT ")"
         " coal_deaths(" SIZE_FORMAT ")"
         " + count(" SSIZE_FORMAT ")",
         p2i(this), size(), _allocation_stats.prev_sweep(), _allocation_stats.split_births(),
         _allocation_stats.coal_births(), _allocation_stats.split_deaths(),
         _allocation_stats.coal_deaths(), count());
}
#endif
template class AdaptiveFreeList<FreeChunk>;
