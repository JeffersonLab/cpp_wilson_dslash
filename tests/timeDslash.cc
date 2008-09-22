#include "unittest.h"
#include "timeDslash.h"

#include "qdp.h"
using namespace QDP;

#undef PAT
#ifdef PAT
#include <pat_api.h>
#endif

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "cpp_dslash.h"
#include "cpp_dslash_qdp_packer.h"

using namespace Assertions;
using namespace std;
using namespace CPlusPlusWilsonDslash;

#ifdef DSLASH_USE_OMP_THREADS
#include <omp.h>
#endif

void
timeDslash::run(void) 
{

  // If we have openmp then do this
#ifdef DSLASH_USE_OMP_THREADS
  int threads_num;
  int myId;

#pragma omp parallel private(threads_num, myId) default(none)
  {
    threads_num = omp_get_num_threads();
    myId = omp_get_thread_num();
    if ( myId == 0 ) { 
      printf("\nRunning with %d OpenMP threads\n", threads_num);
    }
  }
#endif

  LatticeFermionF chi, psi;
  LatticeFermionD chid, psid;

  // Make a random gauge field 
  multi1d<LatticeColorMatrixF> u(4);
  multi1d<LatticeColorMatrixD> ud(4);

  for(int mu=0; mu < 4; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
    ud[mu]=u[mu];
  }

  // Make a random source
  gaussian(psi);
  gaussian(psid);
  
  // Initialize the wilson dslash
  Dslash<float> D32(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);
  // Initialize the wilson dslash
  Dslash<double> D64(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);
  
  /// Pack the gauge fields
  multi1d<PrimitiveSU3MatrixF> packed_gauge;
  packed_gauge.resize( 4 * Layout::sitesOnNode() );
  qdp_pack_gauge(u, packed_gauge);

  multi1d<PrimitiveSU3MatrixD> packed_gauged;
  packed_gauged.resize( 4 * Layout::sitesOnNode() );
  qdp_pack_gauge(ud, packed_gauged);
 

  QDPIO::cout << endl;

  StopWatch swatch;
  double time=0;
  double n_secs = 25;
  int iters=130000;
  QDPIO::cout << endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << endl;


  swatch.reset();
  swatch.start();

#ifdef PAT
  int ierr;
  ierr=PAT_region_begin(19, "DslashLoop");
#endif
  for(int i=0; i < iters; ++i) {
    D32( (float *)&(chi.elem(0).elem(0).elem(0).real()),
	 (float *)&(psi.elem(0).elem(0).elem(0).real()),
	 (float *)&(packed_gauge[0]),
	 1, 0);

  }
#ifdef PAT
  ierr=PAT_region_end(19);
#endif

  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  Internal::globalSum(time);
  time /= (double)Layout::numNodes();

  QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << endl;
  QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << endl;    
  double Mflops = 1320.0f*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
  double perf = Mflops/time;
  QDPIO::cout << "\t Performance is: " << perf << " Mflops (sp) in Total" << endl;
  QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << endl;
  QDPIO::cout << endl;


  QDPIO::cout << "\t Timing with " << iters << " counts" << endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    D64(  (double *)&(chid.elem(0).elem(0).elem(0).real()),
	  (double *)&(psid.elem(0).elem(0).elem(0).real()),
	  (double *)&(packed_gauged[0]),
	   -1, 0);

  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  Internal::globalSum(time);
  time /= (double)Layout::numNodes();

  QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << endl;
  QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << endl;    
  Mflops = 1320.0f*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
  perf = Mflops/time;
  QDPIO::cout << "\t Performance is: " << perf << " Mflops (dp) in Total" << endl;
  QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << endl;

  // Finalize the Dslash
}
