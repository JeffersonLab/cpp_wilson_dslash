#ifndef CPP_DSLASH_MATVEC_64BIT_H
#define CPP_DSLASH_MATVEC_64BIT_H


#include <cpp_dslash_types.h>
#warning "Using C stuff for 64 bit"


namespace CPlusPlusWilsonDslash {
  using namespace Dslash64BitTypes;

    // Gauge Matrix is 3x3x2 (Natural ordering)
    // HalfSpinor is 2x3x2 (spin, color, reim)
    inline 
    void su3_mult(HalfSpinor res, const GaugeMatrix u, const HalfSpinor src) 
    {

      int re=0;
      int im=1;
      for(int spin=0; spin < 2; spin++) { 
	for(int row=0; row < 3; row++) { 
	  res[spin][row][re]=0;
	  res[spin][row][im]=0;
	}

	for(int row=0; row < 3; row++) { 
	  for(int col=0; col < 3; col++) {
	    res[spin][row][re] += u[col][row][re]*src[spin][col][re]
	      -u[col][row][im]*src[spin][col][im];

	    res[spin][row][im] += u[col][row][re]*src[spin][col][im]
	      +u[col][row][im]*src[spin][col][re];


	  }
	}
      }
    }

    inline 
      void su3_adj_mult(HalfSpinor res, GaugeMatrix u, HalfSpinor src)
    {
      int re=0;
      int im=1;

      for(int spin=0; spin < 2; spin++) {
	for(int row=0; row < 3; row++) { 
	  res[spin][row][re] = 0;
	  res[spin][row][im] = 0;
	}
	
	for(int row=0; row < 3; row++) { 
	  for(int col=0; col < 3; col++) {
	    res[spin][row][re] += u[row][col][re]*src[spin][col][re]
	      +  u[row][col][im]*src[spin][col][im];
	    res[spin][row][im] += u[row][col][re]*src[spin][col][im]
	      - u[row][col][im]*src[spin][col][re];
	  }
	}

      }
    }
}

#endif

    



