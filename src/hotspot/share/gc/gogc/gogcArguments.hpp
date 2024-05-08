

#ifndef SHARE_GC_CMS_CMSARGUMENTS_HPP
#define SHARE_GC_CMS_CMSARGUMENTS_HPP

#include "gc/shared/gcArguments.hpp"
#include "gc/shared/genArguments.hpp"

class CollectedHeap;

class GoGCArguments : public GenArguments {
private:
  void disable_adaptive_size_policy(const char* collector_name);
  void set_parnew_gc_flags();

  virtual void initialize();
  virtual CollectedHeap* create_heap();
};

#endif
