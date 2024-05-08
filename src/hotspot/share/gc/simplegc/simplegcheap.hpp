/*
 * simplegcheap.hpp
 *
 *  Created on: 23-Apr-2019
 *      Author: unmesh
 */

#ifndef SRC_HOTSPOT_SHARE_GC_SIMPLEGC_SIMPLEGCHEAP_HPP_
#define SRC_HOTSPOT_SHARE_GC_SIMPLEGC_SIMPLEGCHEAP_HPP_

#include "gc/shared/collectedHeap.hpp"
#include "gc/shared/markBitMap.hpp"
#include "gc/shared/softRefPolicy.hpp"
#include "gc/shared/space.hpp"
#include "services/memoryManager.hpp"
#include "gc/simplegc/simpleGcMonitoringSupport.hpp"
#include "gc/epsilon/epsilonBarrierSet.hpp"
#include "gc/epsilon/epsilon_globals.hpp"
#include "utilities/pair.hpp"

class SimpleGCHeap: public CollectedHeap {
	friend class VMStructs;
private:
	SoftRefPolicy _soft_ref_policy;
	SimpleGCMonitoringSupport* _monitoring_support;
	MemoryPool* _pool;
	GCMemoryManager _memory_manager;
	ContiguousSpace* _space;
	VirtualSpace _virtual_space;
	size_t _max_tlab_size;
	size_t _step_counter_update;
	size_t _step_heap_print;
	int64_t _decay_time_ns;
	volatile size_t _last_counter_update;
	volatile size_t _last_heap_print;
	MemRegion _bitmap_region;
	MarkBitMap _bitmap;

public:
	static SimpleGCHeap* heap();

	SimpleGCHeap() :_memory_manager("SimpleGC Heap"), _space(nullptr) {
	}
	;

	virtual Name kind() const {
		return CollectedHeap::SimpleGC;
	}

	virtual const char* name() const {
		return "SimpleGC";
	}

	virtual SoftRefPolicy* soft_ref_policy() {
		return &_soft_ref_policy;
	}

	virtual jint initialize();
	virtual void post_initialize();
	virtual void initialize_serviceability();

	virtual GrowableArray<GCMemoryManager*> memory_managers();
	virtual GrowableArray<MemoryPool*> memory_pools();

	virtual size_t max_capacity() const {
		return _virtual_space.reserved_size();
	}
	virtual size_t capacity() const {
		return _virtual_space.committed_size();
	}
	virtual size_t used() const {
		return _space->used();
	}

	virtual bool is_in(const void* p) const {
		return _space->is_in(p);
	}

	virtual bool is_scavengable(oop obj) {
		return true;
	}

	virtual bool is_maximal_no_gc() const {
		// No GC is going to happen. Return "we are at max", when we are about to fail.
		return used() == capacity();
	}

	// Allocation
	HeapWord* allocate_work(size_t size);
	HeapWord* allocate_or_collect_work(size_t size);
	virtual HeapWord* mem_allocate(size_t size,
			bool* gc_overhead_limit_was_exceeded);
	virtual HeapWord* allocate_new_tlab(size_t min_size, size_t requested_size,
			size_t* actual_size);

	// TLAB allocation
	virtual bool supports_tlab_allocation() const {
		return true;
	}
	virtual size_t tlab_capacity(Thread* thr) const {
		return capacity();
	}
	virtual size_t tlab_used(Thread* thr) const {
		return used();
	}
	virtual size_t max_tlab_size() const {
		return _max_tlab_size;
	}
	virtual size_t unsafe_max_tlab_alloc(Thread* thr) const;

	virtual void collect(GCCause::Cause cause);
	virtual void do_full_collection(bool clear_all_soft_refs);

	// Heap walking support
	virtual void safe_object_iterate(ObjectClosure* cl);
	virtual void object_iterate(ObjectClosure* cl) {
		safe_object_iterate(cl);
	}


	virtual bool supports_object_pinning() const {
		return false;
	}
	virtual void pin_object(JavaThread* thread, oop obj) {
		
	}
	virtual void unpin_object(JavaThread* thread, oop obj) {
	}

	bool requires_barriers(stackChunkOop obj) const override { return false; }

	bool print_location(outputStream* st, void* addr) const override {return false; }

	// No support for block parsing.
	virtual HeapWord* block_start(const void* addr) const {
		return NULL;
	}
	virtual size_t block_size(const HeapWord* addr) const {
		return 0;
	}
	virtual bool block_is_obj(const HeapWord* addr) const {
		return false;
	}

	// No GC threads
	virtual void print_gc_threads_on(outputStream* st) const {
	}
	virtual void gc_threads_do(ThreadClosure* tc) const {
	}

	// No heap verification
	virtual void prepare_for_verify() {
	}
	virtual void verify(VerifyOption option) {
	}

	virtual jlong millis_since_last_gc() {
		// Report time since the VM start
		return os::elapsed_counter() / NANOSECS_PER_MILLISEC;
	}

	virtual void print_on(outputStream* st) const;
	virtual void print_tracing_info() const;

	void entry_collect(GCCause::Cause cause);

	 // No nmethod handling
	  virtual void register_nmethod(nmethod* nm) {}
	  virtual void unregister_nmethod(nmethod* nm) {}
	  virtual void flush_nmethod(nmethod* nm) {}
	  virtual void verify_nmethod(nmethod* nm) {}

private:
	void print_heap_info(size_t used) const;
	void print_metaspace_info() const;

	void vmentry_collect(GCCause::Cause cause);

	void do_roots(OopClosure* cl, bool everything);
	void process_roots(OopClosure* cl) {
		do_roots(cl, false);
	}
	void process_all_roots(OopClosure* cl) {
		do_roots(cl, true);
	}
	void walk_bitmap(ObjectClosure* cl);
	Pair<char*, size_t> allocate_heap_memory(size_t heap_size_in_bytes, size_t align);
	void update_counters();
	void print_occupancy();
	void init_monitoring_support();
	void allocate_marking_bitmap(Pair<char*, size_t> heap_base_address_and_size);
	void registerBarrier();
	void log_initial_heap_info();
	void init_heap_log_steps(size_t heap_max_byte_size);
	void init_tlab_size();
};

#endif /* SRC_HOTSPOT_SHARE_GC_SIMPLEGC_SIMPLEGCHEAP_HPP_ */
