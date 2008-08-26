#ifndef DISPATCH_SCALAR_H
#define DISPATCH_SCALAR_H

#include <cstdlib>         /* for size_t */
using namespace std;

#include <shift_table_scalar.h>
namespace CPlusPlusWilsonDslash {


  /* This is absolutely horrible -- There must be a better way */

  /* Thread worker argument structure */
  struct ThreadWorkerArgs {
    void *res;           /*!< Spinor either read */
    void *psi;           /*!< Spinor to  write */
    void *u;        /*!< Gauge field - suitably packed */
    void *s;
    int cb;            /*!< Checkerboard (source) */

  };

  /* Functions: Thread Dispatch */
  void dispatchToThreads(void (*func)(size_t, size_t, int, const void *),
			 void* source,
			 void* result, 
			 void* u,
			 void* s,
			 int cb,
			 int n_sites);

}; // namespace

#endif
