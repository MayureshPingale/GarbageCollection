
#include "precompiled.hpp"
#include "gc/cms/GoGCArguments.hpp"
#include "gc/cms/cmsHeap.hpp"
#include "gc/cms/compactibleFreeListSpace.hpp"
#include "gc/shared/cardTableRS.hpp"
#include "gc/shared/gcArguments.hpp"
#include "gc/shared/genCollectedHeap.hpp"
#include "gc/shared/workerPolicy.hpp"
#include "runtime/arguments.hpp"
#include "runtime/globals.hpp"
#include "runtime/globals_extension.hpp"
#include "utilities/defaultStream.hpp"

void GoGCArguments::set_parnew_gc_flags() {
  assert(!UseSerialGC && !UseParallelOldGC && !UseParallelGC && !UseG1GC,
         "control point invariant");
  assert(UseGoGC, "CMS is expected to be on here");

  if (FLAG_IS_DEFAULT(ParallelGCThreads)) {
    FLAG_SET_DEFAULT(ParallelGCThreads, WorkerPolicy::parallel_worker_threads());
    assert(ParallelGCThreads > 0, "We should always have at least one thread by default");
  } else if (ParallelGCThreads == 0) {
    jio_fprintf(defaultStream::error_stream(),
        "The ParNew GC can not be combined with -XX:ParallelGCThreads=0\n");
    vm_exit(1);
  }


  if (FLAG_IS_DEFAULT(YoungPLABSize)) {
    FLAG_SET_DEFAULT(YoungPLABSize, (intx)1024);
  }
  if (FLAG_IS_DEFAULT(OldPLABSize)) {
    FLAG_SET_DEFAULT(OldPLABSize, (intx)1024);
  }

  // When using compressed oops, we use local overflow stacks,
  // rather than using a global overflow list chained through
  // the klass word of the object's pre-image.
  if (UseCompressedOops && !ParGCUseLocalOverflow) {
    if (!FLAG_IS_DEFAULT(ParGCUseLocalOverflow)) {
      warning("Forcing +ParGCUseLocalOverflow: needed if using compressed references");
    }
    FLAG_SET_DEFAULT(ParGCUseLocalOverflow, true);
  }
  assert(ParGCUseLocalOverflow || !UseCompressedOops, "Error");
}


void GoGCArguments::initialize() {
  GCArguments::initialize();

  assert(!UseSerialGC && !UseParallelOldGC && !UseParallelGC, "Error");
  assert(UseGoGC, "GoGC is expected to be on here");


  if (UseGoGC && FLSVerifyAllHeapReferences) {
    if (VerifyDuringStartup) {
      warning("Heap verification at start-up disabled "
              "(due to current incompatibility with FLSVerifyAllHeapReferences)");
      VerifyDuringStartup = false; // Disable verification at start-up
    }

    if (VerifyBeforeExit) {
      warning("Heap verification at shutdown disabled "
              "(due to current incompatibility with FLSVerifyAllHeapReferences)");
      VerifyBeforeExit = false; // Disable verification at shutdown
    }
  }

  if (!ClassUnloading) {
    FLAG_SET_CMDLINE(CMSClassUnloadingEnabled, false);
  }

  // Set CMS global values
  CompactibleFreeListSpace::set_cms_values();

  // Turn off AdaptiveSizePolicy by default for cms until it is complete.
  disable_adaptive_size_policy("UseGoGC");

  set_parnew_gc_flags();

  size_t max_heap = align_down(MaxHeapSize,
                               CardTableRS::ct_max_alignment_constraint());

  // Now make adjustments for CMS
  intx   tenuring_default = (intx)6;
  size_t young_gen_per_worker = CMSYoungGenPerWorker;

  // Preferred young gen size for "short" pauses:
  // upper bound depends on # of threads and NewRatio.
  const size_t preferred_max_new_size_unaligned =
    MIN2(max_heap/(NewRatio+1), ScaleForWordSize(young_gen_per_worker * ParallelGCThreads));
  size_t preferred_max_new_size =
    align_up(preferred_max_new_size_unaligned, os::vm_page_size());

  // Unless explicitly requested otherwise, size young gen
  // for "short" pauses ~ CMSYoungGenPerWorker*ParallelGCThreads

  // If either MaxNewSize or NewRatio is set on the command line,
  // assume the user is trying to set the size of the young gen.
  if (FLAG_IS_DEFAULT(MaxNewSize) && FLAG_IS_DEFAULT(NewRatio)) {

    // Set MaxNewSize to our calculated preferred_max_new_size unless
    // NewSize was set on the command line and it is larger than
    // preferred_max_new_size.
    if (!FLAG_IS_DEFAULT(NewSize)) {   // NewSize explicitly set at command-line
      FLAG_SET_ERGO(MaxNewSize, MAX2(NewSize, preferred_max_new_size));
    } else {
      FLAG_SET_ERGO(MaxNewSize, preferred_max_new_size);
    }
    log_trace(gc, heap)("CMS ergo set MaxNewSize: " SIZE_FORMAT, MaxNewSize);

    // Code along this path potentially sets NewSize and OldSize
    log_trace(gc, heap)("CMS set min_heap_size: " SIZE_FORMAT " initial_heap_size:  " SIZE_FORMAT " max_heap: " SIZE_FORMAT,
                        MinHeapSize, InitialHeapSize, max_heap);
    size_t min_new = preferred_max_new_size;
    if (FLAG_IS_CMDLINE(NewSize)) {
      min_new = NewSize;
    }
    if (max_heap > min_new && MinHeapSize > min_new) {
      // Unless explicitly requested otherwise, make young gen
      // at least min_new, and at most preferred_max_new_size.
      if (FLAG_IS_DEFAULT(NewSize)) {
        FLAG_SET_ERGO(NewSize, MAX2(NewSize, min_new));
        FLAG_SET_ERGO(NewSize, MIN2(preferred_max_new_size, NewSize));
        log_trace(gc, heap)("CMS ergo set NewSize: " SIZE_FORMAT, NewSize);
      }
      // Unless explicitly requested otherwise, size old gen
      // so it's NewRatio x of NewSize.
      if (FLAG_IS_DEFAULT(OldSize)) {
        if (max_heap > NewSize) {
          FLAG_SET_ERGO(OldSize, MIN2(NewRatio*NewSize, max_heap - NewSize));
          log_trace(gc, heap)("CMS ergo set OldSize: " SIZE_FORMAT, OldSize);
        }
      }
    }
  }
  // Unless explicitly requested otherwise, definitely
  // promote all objects surviving "tenuring_default" scavenges.
  if (FLAG_IS_DEFAULT(MaxTenuringThreshold) &&
      FLAG_IS_DEFAULT(SurvivorRatio)) {
    FLAG_SET_ERGO(MaxTenuringThreshold, tenuring_default);
  }
  // If we decided above (or user explicitly requested)
  // `promote all' (via MaxTenuringThreshold := 0),
  // prefer minuscule survivor spaces so as not to waste
  // space for (non-existent) survivors
  if (FLAG_IS_DEFAULT(SurvivorRatio) && MaxTenuringThreshold == 0) {
    FLAG_SET_ERGO(SurvivorRatio, MAX2((uintx)1024, SurvivorRatio));
  }

  // OldPLABSize is interpreted in CMS as not the size of the PLAB in words,
  // but rather the number of free blocks of a given size that are used when
  // replenishing the local per-worker free list caches.
  if (FLAG_IS_DEFAULT(OldPLABSize)) {
    if (!FLAG_IS_DEFAULT(ResizeOldPLAB) && !ResizeOldPLAB) {
      // OldPLAB sizing manually turned off: Use a larger default setting,
      // unless it was manually specified. This is because a too-low value
      // will slow down scavenges.
      FLAG_SET_ERGO(OldPLABSize, CompactibleFreeListSpaceLAB::_default_static_old_plab_size); // default value before 6631166
    } else {
      FLAG_SET_DEFAULT(OldPLABSize, CompactibleFreeListSpaceLAB::_default_dynamic_old_plab_size); // old CMSParPromoteBlocksToClaim default
    }
  }

  // If either of the static initialization defaults have changed, note this
  // modification.
  if (!FLAG_IS_DEFAULT(OldPLABSize) || !FLAG_IS_DEFAULT(OldPLABWeight)) {
    CompactibleFreeListSpaceLAB::modify_initialization(OldPLABSize, OldPLABWeight);
  }

  log_trace(gc)("MarkStackSize: %uk  MarkStackSizeMax: %uk", (unsigned int) (MarkStackSize / K), (uint) (MarkStackSizeMax / K));
}

void GoGCArguments::disable_adaptive_size_policy(const char* collector_name) {
  if (UseAdaptiveSizePolicy) {
    if (FLAG_IS_CMDLINE(UseAdaptiveSizePolicy)) {
      warning("Disabling UseAdaptiveSizePolicy; it is incompatible with %s.",
              collector_name);
    }
    FLAG_SET_DEFAULT(UseAdaptiveSizePolicy, false);
  }
}

CollectedHeap* GoGCArguments::create_heap() {
  return new CMSHeap();
}
