
#ifndef SHARE_GC_gogc_gogc_GLOBALS_HPP
#define SHARE_GC_gogc_gogc_GLOBALS_HPP

#include "runtime/globals_shared.hpp"

//
// Defines all globals flags used by the gogc GC.
//

#define GC_gogc_FLAGS(develop,                                           \
                         develop_pd,                                        \
                         product,                                           \
                         product_pd,                                        \
                         notproduct,                                        \
                         range,                                             \
                         constraint)                                        \
                                                                            \
  product(size_t, gogcPrintHeapSteps, 20, EXPERIMENTAL,                  \
          "Print heap occupancy stats with this number of steps. "          \
          "0 turns the printing off.")                                      \
          range(0, max_intx)                                                \
                                                                            \
  product(size_t, gogcUpdateCountersStep, 1 * M, EXPERIMENTAL,           \
          "Update heap occupancy counters after allocating this much "      \
          "memory. Higher values would make allocations faster at "         \
          "the expense of lower resolution in heap counters.")              \
          range(1, max_intx)                                                \
                                                                            \
  product(size_t, gogcMaxTLABSize, 4 * M, EXPERIMENTAL,                  \
          "Max TLAB size to use with gogc GC. Larger value improves "    \
          "performance at the expense of per-thread memory waste. This "    \
          "asks TLAB machinery to cap TLAB sizes at this value.")           \
          range(1, max_intx)                                                \
                                                                            \
  product(bool, gogcElasticTLAB, true, EXPERIMENTAL,                     \
          "Use elastic policy to manage TLAB sizes. This conserves memory " \
          "for non-actively allocating threads, even when they request "    \
          "large TLABs for themselves. Active threads would experience "    \
          "smaller TLABs until policy catches up.")                         \
                                                                            \
  product(bool, gogcElasticTLABDecay, true, EXPERIMENTAL,                \
          "Use timed decays to shrink TLAB sizes. This conserves memory "   \
          "for the threads that allocate in bursts of different sizes, "    \
          "for example the small/rare allocations coming after the initial "\
          "large burst.")                                                   \
                                                                            \
  product(double, gogcTLABElasticity, 1.1, EXPERIMENTAL,                 \
          "Multiplier to use when deciding on next TLAB size. Larger value "\
          "improves performance at the expense of per-thread memory waste. "\
          "Lower value improves memory footprint, but penalizes actively "  \
          "allocating threads.")                                            \
          range(1.0, DBL_MAX)                                               \
                                                                            \
  product(size_t, gogcTLABDecayTime, 1000, EXPERIMENTAL,                 \
          "TLAB sizing policy decays to initial size after thread had not " \
          "allocated for this long. Time is in milliseconds. Lower value "  \
          "improves memory footprint, but penalizes actively allocating "   \
          "threads.")                                                       \
          range(1, max_intx)                                                \
                                                                            \
  product(size_t, gogcMinHeapExpand, 128 * M, EXPERIMENTAL,              \
          "Min expansion step for heap. Larger value improves performance " \
          "at the potential expense of memory waste.")                      \
          range(1, max_intx)                                                \
                                                                            \
  product(bool, gogcSlidingGC, false, EXPERIMENTAL,                      \
          "Actually does sliding mark-compact GC.")                         \
                                                                            \
  product(bool, gogcImplicitGC, true, EXPERIMENTAL,                      \
          "Does GC on implicit GC requests, e.g. for allocation failure.")  \
                                                                            \
  product(bool, gogcUncommit, false, EXPERIMENTAL,                       \
          "Uncommits all unneeded memory after GC.")                        \
                                                                            \
  product(bool, gogcVerify, false, EXPERIMENTAL,                         \
          "Does the additional GC verification step.")                      \
                                                                            \

// end of GC_gogc_FLAGS

#endif // SHARE_GC_gogc_gogc_GLOBALS_HPP
