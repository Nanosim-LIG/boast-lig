#ifndef __SKB_UTILS__H_
#define __SKB_UTILS__H_

#ifdef SSE
#include "usimd-sse.h"
#elif defined AVX || defined AVX2 
#include "usimd-avx.h"
#elif defined AVX3 || defined LRB1 || defined LRB2
#include "usimd-uisa.h"
#else
#include "usimd-scalar.h"
#endif

#endif
