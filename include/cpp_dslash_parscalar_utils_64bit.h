#ifndef CPP_DSLASH_PARSCALAR_UTILS_64BIT_H
#define CPP_DSLASH_PARSCALAR_UTILS_64BIT_H

#include "dslash_config.h"

#ifdef DSLASH_USE_SSE2
/* SSE DECOMP/RECONS functions */
#include "cpp_dslash_parscalar_decomp_64bit_sse2.h"
#include "cpp_dslash_parscalar_decomp_hvv_64bit_sse2.h"
#include "cpp_dslash_parscalar_mvv_recons_64bit_sse2.h"
#include "cpp_dslash_parscalar_recons_64bit_sse2.h"
#else
/* C Equiv functions */
#include "cpp_dslash_parscalar_decomp_64bit_c.h"
#include "cpp_dslash_parscalar_decomp_hvv_64bit_c.h"
#include "cpp_dslash_parscalar_mvv_recons_64bit_c.h"
#include  "cpp_dslash_parscalar_recons_64bit_c.h"
#endif

#endif
