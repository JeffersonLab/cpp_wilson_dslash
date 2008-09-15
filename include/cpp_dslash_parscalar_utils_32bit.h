#ifndef CPP_DSLASH_PARSCALAR_UTILS_32BIT_H
#define CPP_DSLASH_PARSCALAR_UTILS_32BIT_H

#include "dslash_config.h"


#ifdef DSLASH_USE_SSE2
/* SSE DECOMP/RECONS functions */
#include "cpp_dslash_parscalar_decomp_32bit_sse2.h"
#include "cpp_dslash_parscalar_decomp_hvv_32bit_sse2.h"
#include "cpp_dslash_parscalar_mvv_recons_32bit_sse2.h"
#include "cpp_dslash_parscalar_recons_32bit_sse2.h"
#else 
/* C ones */
#include "cpp_dslash_parscalar_decomp_32bit_c.h"
#include "cpp_dslash_parscalar_decomp_hvv_32bit_c.h"
#include "cpp_dslash_parscalar_mvv_recons_32bit_c.h"
#include "cpp_dslash_parscalar_recons_32bit_c.h"

#endif

#endif
