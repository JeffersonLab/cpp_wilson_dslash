#include <dispatch_parscalar.h>
#include <omp.h>


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
namespace CPlusPlusClover {

void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
		       void* the_spinor,
		       void* the_halfspinor,
		       void *u,
		       void *clov,
		       void *spinor2,
		       void *s,
		       int half,
		       int n_sites)
{
  struct CloverThreadWorkerArgs a;
  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.clov = clov;
  a.spinor2 = spinor2;
  a.s = s;
  a.half = half;


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
}

 
  
}; 
