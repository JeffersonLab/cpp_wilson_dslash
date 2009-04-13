#include "unittest.h"
#include "timeDslash.h"
#include "testvol.h"
#include "cache.h"

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


void
timeDslash::run(void) 
{

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
#if 0
  qdp_pack_gauge(u, packed_gauge);
#else
  // Hand pack the gauge field...
  int volume = Layout::sitesOnNode();
  
  for(int ix = 0; ix < volume; ix++) {
    for(int mu = 0; mu < 4; mu++) {
      packed_gauge[ mu + 4*(ix) ] =
	transpose( u[mu].elem(D32.getPathSite(ix) ).elem() );
    }
  }
#endif


  multi1d<PrimitiveSU3MatrixD> packed_gauged;
  packed_gauged.resize( 4 * Layout::sitesOnNode() );

#if 0
  qdp_pack_gauge(u, packed_gauged);
#else
  // Hand pack the gauge field...

  
  for(int ix = 0; ix < volume; ix++) {
    for(int mu = 0; mu < 4; mu++) {
      packed_gauged[ mu + 4*(ix) ] =
	transpose( ud[mu].elem(D64.getPathSite(ix) ).elem() );
    }
  }
#endif

 

  QDPIO::cout << endl;

  StopWatch swatch;
  double time=0;
  double n_secs = 25;

  QDPIO::cout << endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << endl;




  PrimitiveSpinorF* xpacked_spinor_in;
  PrimitiveSpinorF* xpacked_spinor_out;
  xpacked_spinor_in = (PrimitiveSpinorF*)malloc(volume*sizeof(PrimitiveSpinorF)+Cache::CacheLineSize);
  xpacked_spinor_out = (PrimitiveSpinorF*)malloc(volume*sizeof(PrimitiveSpinorF)+Cache::CacheLineSize);

  if( xpacked_spinor_in == (PrimitiveSpinorF *)NULL ) { 
    cerr << "Fie upon you!" << endl;
    QDP_abort(1);
  }
  if( xpacked_spinor_out == (PrimitiveSpinorF *)NULL ) { 
    cerr << "Fie upon you!" << endl;
    QDP_abort(1);
  }
  ptrdiff_t pad = 0;
  if ( (ptrdiff_t)xpacked_spinor_in % Cache::CacheLineSize != 0 ) {
	pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xpacked_spinor_in % Cache::CacheLineSize);
  }
  PrimitiveSpinorF* packed_spinor_in = (PrimitiveSpinorF *)((char *)xpacked_spinor_in + pad);
 
  pad = 0;
  if ( (ptrdiff_t)xpacked_spinor_out % Cache::CacheLineSize != 0 ) {
	pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xpacked_spinor_out % Cache::CacheLineSize);
  }
  PrimitiveSpinorF* packed_spinor_out = (PrimitiveSpinorF *)((char *)xpacked_spinor_out + pad);
  


  for(int ix=0; ix < volume; ix++) { 
    packed_spinor_in[ ix ] = psi.elem( D32.getPathSite(ix) ) ;
  }

  swatch.reset();
  swatch.start();
#ifdef PAT
  int ierr;
  ierr=PAT_region_begin(19, "DslashLoop");
#endif


  for(int i=0; i < iters; ++i) {
    D32( (float *)&(packed_spinor_out[0].elem(0).elem(0).real()),
	 (float *)&(packed_spinor_in[0].elem(0).elem(0).real()),
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
  
  free(xpacked_spinor_in);
  free(xpacked_spinor_out);

  QDPIO::cout << "\t Timing with " << iters << " counts" << endl;

  PrimitiveSpinorD* xpacked_spinor_in_d;
  PrimitiveSpinorD* xpacked_spinor_out_d;
  xpacked_spinor_in_d = (PrimitiveSpinorD*)malloc(volume*sizeof(PrimitiveSpinorD)+Cache::CacheLineSize);
  xpacked_spinor_out_d = (PrimitiveSpinorD*)malloc(volume*sizeof(PrimitiveSpinorD)+Cache::CacheLineSize);

  if( xpacked_spinor_in_d == (PrimitiveSpinorD *)NULL ) { 
    cerr << "Fie upon you!" << endl;
    QDP_abort(1);
  }
  if( xpacked_spinor_out_d == (PrimitiveSpinorD *)NULL ) { 
    cerr << "Fie upon you!" << endl;
    QDP_abort(1);
  }
  pad = 0;
  if ( (ptrdiff_t)xpacked_spinor_in_d % Cache::CacheLineSize != 0 ) {
	pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xpacked_spinor_in_d % Cache::CacheLineSize);
  }
  PrimitiveSpinorD* packed_spinor_in_d = (PrimitiveSpinorD *)((char *)xpacked_spinor_in_d + pad);
 
  pad = 0;
  if ( (ptrdiff_t)xpacked_spinor_out_d % Cache::CacheLineSize != 0 ) {
	pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xpacked_spinor_out_d % Cache::CacheLineSize);
  }
  PrimitiveSpinorD* packed_spinor_out_d = (PrimitiveSpinorD *)((char *)xpacked_spinor_out_d + pad);
  


  for(int ix=0; ix < volume; ix++) { 
    packed_spinor_in_d[ ix ] = psid.elem( D64.getPathSite(ix) ) ;
  }




  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    D64(  (double *)&(packed_spinor_out_d[0].elem(0).elem(0).real()),
	  (double *)&(packed_spinor_in_d[0].elem(0).elem(0).real()),
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
  free(xpacked_spinor_in_d);
  free(xpacked_spinor_out_d);

}
