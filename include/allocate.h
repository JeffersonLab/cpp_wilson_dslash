/*
 * allocate.h
 *
 *  Created on: Aug 26, 2016
 *      Author: bjoo
 */

#ifndef CPP_WILSON_DSLASH_INCLUDE_ALLOCATE_H_
#define CPP_WILSON_DSLASH_INCLUDE_ALLOCATE_H_

#include "dslash_config.h"

#ifdef SSE_USE_QDPXX
#include "qdp_pool_allocator.h"
#else
#include <cstdlib>
#endif

namespace  CPlusPlusWilsonDslash {

inline
void *alloc(size_t numbytes) {
#ifdef SSE_USE_QDPXX
	return (void *)QDP::Allocator::theQDPPoolAllocator::Instance().alloc(numbytes);
#else
	return (void *)std::malloc(numbytes);
#endif
}

inline
void dealloc(void *mem) {
#ifdef SSE_USE_QDPXX
	QDP::Allocator::theQDPPoolAllocator::Instance().free(mem);
#else
	std::free(mem);
#endif
}

}


#endif /* OTHER_LIBS_CPP_WILSON_DSLASH_INCLUDE_ALLOCATE_H_ */
