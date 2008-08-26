#include <dispatch_parscalar.h>
#include <qmt.h>

namespace CPlusPlusWilsonDslash {

  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* the_spinor,
			 void* the_halfspinor, 
			 void *u,
			 void *s,
			 int cb,
			 int n_sites)
  {
    ThreadWorkerArgs a;
    int threads_num;
    int chucksize;
    int myId;
    int low;
    int high;
    
    a.spinor = the_spinor;
    a.half_spinor = the_halfspinor;
    a.u = u;
    a.cb = cb;
    a.s = s;
    qmt_call((qmt_userfunc_t)func, n_sites, &a);
  }

}; // End Namespace
