#include <dispatch_scalar.h>
#include <omp.h>


namespace CPlusPlusWilsonDslash {
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* source,
			 void* result, 
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
   
    a.psi = source;
    a.res = result;
    a.u = u;
    a.cb = cb;
    a.s = s;
#pragma omp parallel shared(func, n_sites, a)				\
  private(threads_num, chucksize, myId, low, high) default(none)
    {
      
      threads_num = omp_get_num_threads();
      chucksize = n_sites/threads_num;
      myId = omp_get_thread_num();
      low = chucksize * myId;
      high = chucksize * (myId+1);
      (*func)(low, high, myId, &a);
    }
  }

}; // End Namespace


// Clover dispatcher
namespace CPlusPlusClover { 

  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* source,
			 void* result, 
			 void* u,
			 void* invclov_ee,
			 void* clov_oo,
			 void* t_spinor,
			 void* s,
			 int n_sites) 
  {
  struct CloverThreadWorkerArgs a;
  a.psi = source;
  a.res =result;
  a.u = u;
  a.invclov_ee = invclov_ee;
  a.clov_oo = clov_oo;
  a.t_spinor = t_spinor;
  a.s = s;


  int threads_num;
  int chucksize;
  int myId;
  int low;
  int high;
#pragma omp parallel shared(func, n_sites, a)				\
  private(threads_num, chucksize, myId, low, high) default(none)
    {
      
      threads_num = omp_get_num_threads();
      chucksize = n_sites/threads_num;
      myId = omp_get_thread_num();
      low = chucksize * myId;
      high = chucksize * (myId+1);
      (*func)(low, high, myId, &a);
    }
  } 

} // End namespace
