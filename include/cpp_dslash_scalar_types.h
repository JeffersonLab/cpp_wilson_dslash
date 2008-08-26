#ifndef CPP_DSLASH_SCALAR_TYPES_H
#define CPP_DSLASH_SCALAR_TYPES_H

namespace CPlusPlusWilsonDslash {



  namespace DslashScalar32BitTypes {
   typedef float FourSpinor[4][3][2];
   typedef float HalfSpinor[3][2][2];
   typedef float GaugeMatrix[3][3][2];
  }

  namespace DslashScalar64BitTypes {
   typedef double FourSpinor[4][3][2];
   typedef double HalfSpinor[2][3][2];
   typedef double GaugeMatrix[3][3][2];
  }


}


#endif
