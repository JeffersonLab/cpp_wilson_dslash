#ifndef CPP_DSLASH_QDP_PACKER_H
#define CPP_DSLASH_QDP_PACKER_H

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif 

namespace CPlusPlusWilsonDslash { 

  typedef PColorMatrix<RComplex<REAL32>, 3> PrimitiveSU3MatrixF;
  typedef PColorMatrix<RComplex<REAL64>, 3> PrimitiveSU3MatrixD;


  void qdp_pack_gauge(const multi1d<LatticeColorMatrixF>&_u, multi1d<PrimitiveSU3MatrixF>& u_tmp);
  
  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrixF>&_u, multi1d<PrimitiveSU3MatrixF>& u_tmp);

  void qdp_pack_gauge(const multi1d<LatticeColorMatrixD>&_u, multi1d<PrimitiveSU3MatrixD>& u_tmp);

  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrixD>&_u, multi1d<PrimitiveSU3MatrixD>& u_tmp);

};

#endif
