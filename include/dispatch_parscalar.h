#ifndef DISPATCH_PARSCALAR_H
#define DISPATCH_PARSCALAR_H

#include <cstdlib>         /* for size_t */
using namespace std;

namespace CPlusPlusWilsonDslash {


  /* This is absolutely horrible -- There must be a better way */

  /* Thread worker argument structure */
  struct ThreadWorkerArgs {
    void *spinor;           /*!< Spinor either read */
    void *half_spinor;           /*!< Spinor to  write */
    void *u;        /*!< Gauge field - suitably packed */
    void *s;
    int cb;            /*!< Checkerboard (source) */
  };

  /* Functions: Thread Dispatch */
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* the_spinor,
			 void* the_halfspinor, 
			 void* u,
			 void* s,
			 int cb,
			 int n_sites);

}; // namespace

#endif
