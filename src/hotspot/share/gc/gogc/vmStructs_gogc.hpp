#ifndef SHARE_GC_gogc_VMSTRUCTS_gogc_HPP
#define SHARE_GC_gogc_VMSTRUCTS_gogc_HPP

#include "gc/gogc/gogcHeap.hpp"
#include "gc/gogc/compactibleFreeListSpace.hpp"
#include "gc/gogc/concurrentMarkSweepGeneration.hpp"
#include "gc/gogc/concurrentMarkSweepThread.hpp"
#include "gc/gogc/parNewGeneration.hpp"

#define VM_STRUCTS_gogcGC(nonstatic_field,                                                                                            \
                         volatile_nonstatic_field,                                                                                   \
                         static_field)                                                                                               \
  nonstatic_field(CompactibleFreeListSpace,    _collector,                                    gogcCollector*)                         \
  nonstatic_field(CompactibleFreeListSpace,    _bt,                                           BlockOffsetArrayNonContigSpace)        \
     static_field(CompactibleFreeListSpace,    _min_chunk_size_in_bytes,                      size_t)                                \
  nonstatic_field(gogcBitMap,                   _bmStartWord,                                  HeapWord*)                             \
  nonstatic_field(gogcBitMap,                   _bmWordSize,                                   size_t)                                \
  nonstatic_field(gogcBitMap,                   _shifter,                                      const int)                             \
  nonstatic_field(gogcBitMap,                   _bm,                                           BitMapView)                            \
  nonstatic_field(gogcBitMap,                   _virtual_space,                                VirtualSpace)                          \
  nonstatic_field(gogcCollector,                _markBitMap,                                   gogcBitMap)                             \
  nonstatic_field(ConcurrentMarkSweepGeneration, _gogcSpace,                                   CompactibleFreeListSpace*)             \
     static_field(ConcurrentMarkSweepThread,   _collector,                                    gogcCollector*)                         \
  nonstatic_field(LinearAllocBlock,            _word_size,                                    size_t)                                \
  nonstatic_field(AFLBinaryTreeDictionary,     _total_size,                                   size_t)                                \
  nonstatic_field(CompactibleFreeListSpace,    _dictionary,                                   AFLBinaryTreeDictionary*)              \
  nonstatic_field(CompactibleFreeListSpace,    _indexedFreeList[0],                           AdaptiveFreeList<FreeChunk>)           \
  nonstatic_field(CompactibleFreeListSpace,    _smallLinearAllocBlock,                        LinearAllocBlock)                      \
  volatile_nonstatic_field(FreeChunk,          _size,                                         size_t)                                \
  nonstatic_field(FreeChunk,                   _next,                                         FreeChunk*)                            \
  nonstatic_field(FreeChunk,                   _prev,                                         FreeChunk*)                            \
  nonstatic_field(AdaptiveFreeList<FreeChunk>, _size,                                         size_t)                                \
  nonstatic_field(AdaptiveFreeList<FreeChunk>, _count,                                        ssize_t)



#define VM_TYPES_gogcGC(declare_type,                                      \
                       declare_toplevel_type,                             \
                       declare_integer_type)                              \
                                                                          \
           declare_type(gogcHeap,                      GenCollectedHeap)   \
           declare_type(ConcurrentMarkSweepGeneration,CardGeneration)     \
           declare_type(ParNewGeneration,             DefNewGeneration)   \
           declare_type(CompactibleFreeListSpace,     CompactibleSpace)   \
           declare_type(ConcurrentMarkSweepThread,    NamedThread)        \
  declare_toplevel_type(gogcCollector)                                     \
  declare_toplevel_type(gogcBitMap)                                        \
  declare_toplevel_type(FreeChunk)                                        \
  declare_toplevel_type(metaspace::Metablock)                             \
  declare_toplevel_type(ConcurrentMarkSweepThread*)                       \
  declare_toplevel_type(ConcurrentMarkSweepGeneration*)                   \
  declare_toplevel_type(CompactibleFreeListSpace*)                        \
  declare_toplevel_type(gogcCollector*)                                    \
  declare_toplevel_type(AFLBinaryTreeDictionary)                          \
  declare_toplevel_type(LinearAllocBlock)                                 \
  declare_toplevel_type(FreeChunk*)                                       \
  declare_toplevel_type(AdaptiveFreeList<FreeChunk>*)                     \
  declare_toplevel_type(AdaptiveFreeList<FreeChunk>)


#define VM_INT_CONSTANTS_gogcGC(declare_constant,                          \
                               declare_constant_with_value)               \
  declare_constant(CompactibleFreeListSpace::IndexSetSize)                \
  declare_constant(Generation::ConcurrentMarkSweep)                       \
  declare_constant(Generation::ParNew)

#endif // SHARE_GC_gogc_VMSTRUCTS_gogc_HPP
