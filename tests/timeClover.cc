#include "unittest.h"
#include "timeClover.h"

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

#include "cpp_clover_scalar.h"
#include "cpp_dslash_qdp_packer.h"

using namespace Assertions;
using namespace std;
using namespace CPlusPlusClover;

#ifdef DSLASH_USE_OMP_THREADS
#include <omp.h>
#endif

void
timeClover::run(void) 
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

  LatticeFermionF3 chi, psi;
  LatticeFermionD3 chid, psid;

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

  int iters=32000;

  {

    CloverSchur4D<float> Klov32(Layout::lattSize().slice(),
				Layout::QDPXX_getSiteCoords,
				Layout::QDPXX_getLinearSiteIndex,
				Layout::QDPXX_nodeNumber);
    
    
    
    /// Pack the gauge fields
    multi1d<PrimitiveSU3MatrixF> packed_gauge;
    packed_gauge.resize( 4 * Layout::sitesOnNode() );
    qdp_pack_gauge(u, packed_gauge);
    
    
    multi1d<Clover32BitTypes::CloverTerm> clov __attribute__((aligned(16)));
    clov.resize(Layout::sitesOnNode());
    multi1d<Clover32BitTypes::CloverTerm> invclov __attribute__((aligned(16)));
    invclov.resize(Layout::sitesOnNode());
    
    // Randomize clover term....
    for(int site=0; site < Layout::sitesOnNode(); site++) { 
      for(int j=0; j < 2; j++) { 
	for(int d=0; d < 6; d++) { 
	  (clov[site])[j].diag[d] = drand48();
	  (invclov[site])[j].diag[d] = drand48();
	}
	for(int od=0; od < 15; od++) { 
	  (clov[site])[j].off_diag[od][0] = drand48();
	  (clov[site])[j].off_diag[od][1] = drand48();
	  (invclov[site])[j].off_diag[od][0] = drand48();
	  (invclov[site])[j].off_diag[od][1] = drand48();
	}
      }
    }
    
    QDPIO::cout << endl;
    
    StopWatch swatch;
    QDPIO::cout << "\t Timing with " << iters << " counts" << endl;


    swatch.reset();
    swatch.start();

    for(int i=0; i < iters; i++) { 
      // Apply clover op to psi into chi
      Klov32((float *)&(chi.elem(0).elem(0).elem(0).real()),	  
	     (float *)&(psi.elem(0).elem(0).elem(0).real()),
	     (float *)&(packed_gauge[0]),
	     (float *)&(invclov[0]),
	     (float *)&(clov[0]),
	     1); 
    }
  
    swatch.stop();
    double time=swatch.getTimeInSeconds();
    
    // Average time over nodes
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << endl;
    QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << endl;    
    // Mflops=2Dslash + 2Clover Terms + 24 flops for subtractig them
    double Mflops =2.0f*(552.0f + 1320.0f+24.0f)*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
    double perf = Mflops/time;
    QDPIO::cout << "\t Performance is: " << perf << " Mflops (sp) in Total" << endl;
    QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << endl;
    QDPIO::cout << endl;
  }




  { 
    CloverSchur4D<double> Klov64(Layout::lattSize().slice(),
				 Layout::QDPXX_getSiteCoords,
				 Layout::QDPXX_getLinearSiteIndex,
				 Layout::QDPXX_nodeNumber);
    

    // Pack Gauge field
    multi1d<PrimitiveSU3MatrixD> packed_gauged;
    packed_gauged.resize( 4 * Layout::sitesOnNode() );
    qdp_pack_gauge(ud, packed_gauged);
    
    multi1d<Clover64BitTypes::CloverTerm> clovd __attribute__((aligned(16)));
    clovd.resize(Layout::sitesOnNode());
    multi1d<Clover64BitTypes::CloverTerm> invclovd __attribute__((aligned(16)));
    invclovd.resize(Layout::sitesOnNode());
    
    // Randomize dprec clover term....
    for(int site=0; site < Layout::sitesOnNode(); site++) { 
      for(int j=0; j < 2; j++) { 
	for(int d=0; d < 6; d++) { 
	  (clovd[site])[j].diag[d] = drand48();
	  (invclovd[site])[j].diag[d] = drand48();
	}
	for(int od=0; od < 15; od++) { 
	  (clovd[site])[j].off_diag[od][0] = drand48();
	  (clovd[site])[j].off_diag[od][1] = drand48();
	  (invclovd[site])[j].off_diag[od][0] = drand48();
	  (invclovd[site])[j].off_diag[od][1] = drand48();
	}
      }
    }
    
    QDPIO::cout << "\t Timing with " << iters << " counts" << endl;
    
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    

    for(int i=0; i < iters; i++) { 
      // Apply clover op to psi into chi
      Klov64((double *)&(chid.elem(0).elem(0).elem(0).real()),	  
	     (double *)&(psid.elem(0).elem(0).elem(0).real()),
	     (double *)&(packed_gauged[0]),
	     (double *)&(invclovd[0]),
	     (double *)&(clovd[0]),
	     1); 
    }
  
    swatch.stop();
    double time=swatch.getTimeInSeconds();
    
    // Average time over nodes
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
  
    QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << endl;
    QDPIO::cout << "\t " << 1.0e6*time/(double)iters << " u sec/iteration" << endl;    
    // Mflops=2Dslash + 2Clover Terms + 24 flops for subtractig them
    double Mflops =2.0f*(552.0f + 1320.0f+24.0f)*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
    double perf = Mflops/time;
    QDPIO::cout << "\t Performance is: " << perf << " Mflops (sp) in Total" << endl;
    QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << endl;
    QDPIO::cout << endl;
  }

}
