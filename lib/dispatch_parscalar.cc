#include <dispatch_parscalar.h>

namespace CPlusPlusWilsonDslash {

void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
		       void* the_spinor,
		       void* the_halfspinor,
		       void *u,
		       void *s,
		       int cb,
		       int n_sites)
{
  struct ThreadWorkerArgs a;
  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.cb = cb;
  a.s = s;
  /* Call dispatch function, with lo=0, hi=n_sites */
  (*func)(0, n_sites, 0, &a);
}

  
}; // End Namespace
