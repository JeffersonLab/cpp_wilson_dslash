#include "unittest.h"
#include "testDslashFull.h"
#include "cache.h"

#include "qdp.h"
using namespace QDP;

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
testDslashFull::run(void) 
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

  LatticeFermionF3 chi, chi2, psi;
  LatticeFermionD3 chid, chi2d, psid;

  // What we consider to be small enough...
  Double small32,small64;
  small32 = Double(1.0e-7);
  small64 = Double(1.0e-16);


  // Make a random gauge field 
  multi1d<LatticeColorMatrixF> u(4);
  multi1d<LatticeColorMatrixD> ud(4);

  for(int mu=0; mu < 4; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
    ud[mu] = u[mu];
  }

  // Make a random source
  gaussian(psi);
  // Downcast
  psid=psi;

  Dslash<float> D32(Layout::lattSize().slice(),
		    Layout::QDPXX_getSiteCoords,
		    Layout::QDPXX_getLinearSiteIndex,
		    Layout::QDPXX_nodeNumber);
  

  //  QDPIO::cout << "Need to allocate the packed gauge: "<< 4*Layout::sitesOnNode()*sizeof(PrimitiveSU3MatrixF) << endl << flush;
  
  multi1d<PrimitiveSU3MatrixF> packed_gauge __attribute__((aligned(16)));
  packed_gauge.resize(4*Layout::sitesOnNode());

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

  // Spinor in
  multi1d<PrimitiveSpinorF> packed_spinor_in  __attribute__((aligned(16)));
  multi1d<PrimitiveSpinorF> packed_spinor_out __attribute__((aligned(16)));
  packed_spinor_in.resize(Layout::sitesOnNode());
  packed_spinor_out.resize(Layout::sitesOnNode());
  for(int ix=0; ix < volume; ix++) { 
    packed_spinor_in[ ix ] = psi.elem( D32.getPathSite(ix) ) ;
  }

#endif
  QDPIO::cout << endl;

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chi = zero;
      for(int ix=0; ix < volume; ix++) { 
	packed_spinor_out[ ix ] = chi.elem( D32.getPathSite(ix) ) ;
      }

      chi2 = zero;

      // Apply SSE Dslash
      D32((float *)&(packed_spinor_out[0].elem(0).elem(0).real()),	
	  (float *)&(packed_spinor_in[0].elem(0).elem(0).real()),
	  (float *)&(packed_gauge[0]),
	  isign, 
	  source_cb);
      
      // Apply QDP Dslash
      dslash(chi2,u,psi, isign, target_cb);
      
      // Export Fermion
      for(int ix=0; ix < volume; ix++) { 
	chi.elem( D32.getPathSite(ix) ) = packed_spinor_out[ ix ] ;
      }


      // Check the difference per number in chi vector
      LatticeFermion diff = chi2 -chi;

      Double diff_norm = sqrt( norm2( diff ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << source_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      assertion( toBool( diff_norm < small32 ) );

    }
  }

  Dslash<double> D64(Layout::lattSize().slice(),
			 Layout::QDPXX_getSiteCoords,
			 Layout::QDPXX_getLinearSiteIndex,
			 Layout::QDPXX_nodeNumber);


   /// Pack the gauge fields
   multi1d<PrimitiveSU3MatrixD> packed_gauged __attribute__((aligned(16)));
   packed_gauged.resize( 4 * Layout::sitesOnNode() );

#if 0
  qdp_pack_gauge(ud, packed_gauged);
#else
  // Hand pack the gauge field...
  volume = Layout::sitesOnNode();
  
  for(int ix = 0; ix < volume; ix++) 
  {
    for(int mu = 0; mu < 4; mu++) 
    { 
      packed_gauged[ mu + 4*(ix) ] =
	transpose( ud[mu].elem(D32.getPathSite(ix)).elem() );
    }
  }
 // Spinor in
  //  multi1d<PrimitiveSpinorD> packed_spinor_in_d  __attribute__((aligned(16)));
  //  multi1d<PrimitiveSpinorD> packed_spinor_out_d __attribute__((aligned(16)));
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
  ptrdiff_t pad = 0;
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

#endif

   // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
   for(int isign=1; isign >= -1; isign -=2) {
     for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chid = zero;
      for(int ix=0; ix < volume; ix++) { 
	packed_spinor_out_d[ ix ] = chid.elem( D64.getPathSite(ix) ) ;
      }

      chi2 = zero;

      // Apply SSE Dslash
      D64((double *)&(packed_spinor_out_d[0].elem(0).elem(0).real()),	  
	  (double *)&(packed_spinor_in_d[0].elem(0).elem(0).real()),
	  (double *)&(packed_gauged[0]),
	  isign, 
	  source_cb);
      
      // Apply QDP Dslash


      dslash(chi2,u,psi, isign, target_cb);

      // Export Fermion
      for(int ix=0; ix < volume; ix++) { 
	chid.elem( D64.getPathSite(ix) ) = packed_spinor_out_d[ ix ] ;
      }


      chi = chid;

      
      // Check the difference per number in chi vector
      LatticeFermionF diff = chi2 - chi;

      Double diff_norm = sqrt( norm2( diff ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << source_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      assertion( toBool( diff_norm < small32 ) );

    }
   }
   free(xpacked_spinor_in_d);
   free(xpacked_spinor_out_d);

}
