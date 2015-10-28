#ifndef __USIMD_256
#define __USIMD_256

#include <stdio.h>


//#include <ia32intrin.h>
#include <immintrin.h>
#ifdef _WIN32
#include <intrin.h>
#endif

#define SIMDTYPE_I32 __m256i // only for AVX2

#ifdef AVX2
  #define _MM_LOADU_SI _mm256_loadu_si256
  #define _MM_LOAD_SI(A) _mm256_load_si256((__m256i *)A)
  #define _MM_STOREU_SI _mm256_storeu_si256
  #define _MM_STORE_SI(A,B) _mm256_store_si256((__m256i *)A,B)
  #define _MM_SETZERO_SI _mm256_setzero_si256
  #define _MM_MUL_SI _mm256_mullo_epi32
  #define _MM_ADD_SI _mm256_add_epi32
  #define _MM_LOAD1_SI(A) _mm256_castps_si256(_mm256_broadcast_ss(A))
  #define _MM_CMPEQ_SI _mm256_cmpeq_epi32
  #define _MM_CMPGT_SI _mm256_cmpgt_epi32
  #define _MM_OR_SI _mm256_or_si256
  #define _MM_CASTINT_FLOAT _mm256_castsi256_ps
#endif

#define SIMD_ALIGN 256



#if defined(PRECISION) && (PRECISION == 1)
  #define fptype float
  #define fpconst(n) (n ## f)  
  #define MPI_FPTYPE MPI_FLOAT
  #define VLEN 8
//  #define NNZ_BITS (double)(4*8)
  #define SIMDFPTYPE __m256
  #define SIMDMASK  SIMDFPTYPE
  #undef CEIL_SIMD
  #define CEIL_SIMD(A)  ((int) ceil((double)(A)/(double)VLEN)*VLEN)

//  #define vRngGaussian vsRngGaussian
//  #define vExp vsExp
//  #define trmm strmm
//  #define vInv vsInv
 
  #define _MM_CASTINT_PS _mm256_castsi256_ps
  #define _MM_SETZERO _mm256_setzero_ps
  #define _MM_SET1 _mm256_set1_ps
  #define _MM_BLENDV _mm256_blendv_ps
  #define _MM_BLEND _mm256_blend_ps
  #define _MM_LOAD _mm256_load_ps
  #define _MM_LOADINT(a) _mm256_load_si256((__m256i*)a)
  #define _MM_LOAD1 _mm256_broadcast_ss
  #define _MM_LOAD1INT(a) _mm256_castps_si256(_mm256_broadcast_ss((float *)a))
  #define _MM_LOADU1arg _MM_LOADU
  #define _MM_LOADU _mm256_loadu_ps
  //__m256 _MM_LOAD4(float *addr) {__m128 x = _mm_load_ps(addr); return _mm256_insertf128_ps(_mm256_castps128_ps256(x), x, 1); }
  #define _MM_STORE _mm256_store_ps
  #define _MM_STORE1 _mm_store_ss
  #define _MM_STOREU _mm256_storeu_ps
 
  #define _MM_RSQRT _mm256_rsqrt_ps
  #define _MM_SQRT _mm256_sqrt_ps
  #define _MM_FLOOR _mm256_floor_ps
  #define _MM_MUL _mm256_mul_ps
  #define _MM_ADD _mm256_add_ps
  #define _MM_SUB _mm256_sub_ps
  #define _MM_MAX _mm256_max_ps
  #define _MM_MIN _mm256_min_ps
  #define _MM_AND _mm256_and_ps
  #define _MM_ANDNOT _mm256_andnot_ps
  #define _MM_OR _mm256_or_ps
  #define _MM_XOR _mm256_xor_ps
  #define _MM_DIV _mm256_div_ps
  #define _MM_RCP _mm256_rcp_ps
  #define _MM_EXP _mm256_exp_ps
  #define _MM_LOG _mm256_log_ps
  #define _MM_CMPEQ(a, b) _mm256_cmp_ps(a, b, _CMP_EQ_UQ)
  #define _MM_CMPNEQ(a, b) _mm256_cmp_ps(a, b, _CMP_NEQ_UQ)
  #define _MM_CMPLT(a, b) _mm256_cmp_ps(a, b, _CMP_LT_OS)
  #define _MM_CMPLE(a, b) _mm256_cmp_ps(a, b, _CMP_LE_OS)
  #define _MM_CMPNLE(a, b) _mm256_cmp_ps(a, b, _CMP_NLE_US)
  #define _MM_CMPNLT(a, b) _mm256_cmp_ps(a, b, _CMP_NLT_US)
  #define _MM_INTCMPGT _mm256_cmpgt_epi32
  #define _MM_INTCMPEQ _mm256_cmpeq_epi32
  #define _MM_INTOR _mm256_or_si256
  #define _MM_HADD _mm256_hadd_ps
  #define _MM_DP _mm256_dp_ps
  #define _MM_BLENDV _mm256_blendv_ps
//  extern "C" { __m256 vmlsExp4(__m256 a); }
  
  inline SIMDFPTYPE _MM_BROADCAST (const float *a) {return _mm256_broadcast_ss(a); }
  inline SIMDFPTYPE _MM_BROADCAST_V (const float a) { return _mm256_broadcast_ss(&a); }

  #define _MM_CAST_TO_SI _mm256_castps_si256
  #define _MM_CAST_FROM_SI _mm256_castsi256_ps

  #define _MM_SPLAT0(xmm_reg) _mm256_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(0,0,0,0))
  #define _MM_SPLAT1(xmm_reg) _mm256_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(1,1,1,1))
  #define _MM_SPLAT2(xmm_reg) _mm256_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(2,2,2,2))
  #define _MM_SPLAT3(xmm_reg) _mm256_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(3,3,3,3))

#ifdef AVX2
  #define _MM_MADD213 			_mm256_fmadd_ps
  #define _MM_MSUB213 			_mm256_fmsub_ps
  #define _MM_MADD132(a, dest, b) 	_mm256_fmadd_ps(a, b, dest)		
  #define _MM_MSUB132(a, dest, b) 	_mm256_fmsub_ps(a, b, dest)		
  #define _MM_MADD231(dest, a, b) 	_mm256_fmadd_ps(a, b, dest)
  #define _MM_MSUB231(dest, a, b) 	_mm256_fmsub_ps(a, b, dest)	// a*b-c
  #define _MM_MSUBR231(a, b, c)	_MM_SUB(a, _MM_MUL(b, c))
  #define _MM_ABS(A)                    _mm256_castsi256_ps(_mm256_abs_epi32(_mm256_castps_si256(A)))
#else
  #define _MM_MADD213(a,b,c)	_MM_ADD(c, _MM_MUL(a, b))
  #define _MM_MSUB213(a,b,c)	_MM_SUB(_MM_MUL(a, b),c)
  #define _MM_MADD132(a, b, c) 	_MM_ADD(b, _MM_MUL(a, c))
  #define _MM_MSUB132(a, b, c) 	_MM_SUB(_MM_MUL(a, c), b)
  #define _MM_MADD231(a, b, c) 	_MM_ADD(a, _MM_MUL(b, c))
  #define _MM_MSUB231(a, b, c) 	_MM_SUB(_MM_MUL(b, c), a)
  #define _MM_MSUBR231(a, b, c)	_MM_SUB(a, _MM_MUL(b, c))
#endif

  #define _MM_MOVEMASK _mm256_movemask_ps
  #define _MM_PERMUTE _mm256_permute_ps
  #define _MM_PERMUTE2F128 _mm256_permute2f128_ps
  #define _MM_INSERTF128 _mm256_insertf128_ps
  #define _MM_EXTRACTF128 _mm256_extractf128_ps
  #define _MM_SHUFPS _mm256_shuffle_ps
  #define _MM_CAST128TO256 _mm256_castps128_ps256
  #define _MM_CAST256TO128 _mm256_castps256_ps128
  #define _MM_CASTSI_PS _mm256_castsi256_ps
  #define _MM_CASTPS_SI _mm256_castps_si256  
  #define _MM_UNPACK_LO _mm256_unpacklo_ps
  #define _MM_UNPACK_HI _mm256_unpackhi_ps

  #define ABS fabsf
  #define _MM_GATHER _mm256_gather_ps_nomask
  #define _MM_GATHER_Scalar(base, ind, dense, vdense) dense[0]=base[ind[0]]; dense[1]=base[ind[1]];\
                                                      dense[2]=base[ind[2]]; dense[3]=base[ind[3]]; \
                                                      dense[4]=base[ind[4]]; dense[5]=base[ind[5]]; \
                                                      dense[6]=base[ind[6]]; dense[7]=base[ind[7]]; \
                                                      vdense = _MM_LOAD(dense);
                                        

  #define MONE -1.0f
  #define ZERO 0.0f
  #define ONE 1.0f
  #define TWO 2.0f
  #define THRESH 0.0001f

  __declspec(align(64)) static int MASK[9][8]={
                         {0,0,0,0,0,0,0,0}, 
                         {0x80000000,0,0,0,0,0,0,0}, 
                         {0x80000000,0x80000000,0,0,0,0,0,0},
                         {0x80000000,0x80000000,0x80000000,0,0,0,0,0}, 
                         {0x80000000,0x80000000,0x80000000,0x80000000,0,0,0,0},
                         {0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0,0,0},
                         {0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0,0}, 
                         {0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0}, 
                         {0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000}};


#define _mm_vprefetch1(a, b) ;
#define _mm_vprefetch2(a, b) ;
#define fast_gather_ps(base,idx,output) \
	{ \
		__m128 x1, x2; \
		unsigned long long index; \
		unsigned idx0,idx2,idx4,idx6; \
		unsigned long long idx1,idx3,idx5,idx7; \
	\
		index = ((unsigned long long *)(idx))[0]; \
		idx0 = index; \
		idx1 = index>>32; \
		index = ((unsigned long long *)(idx))[1]; \
		idx2 = index; \
		idx3 = index>>32; \
		index = ((unsigned long long *)(idx))[2]; \
		idx4 = index; \
		idx5 = index>>32; \
		index = ((unsigned long long *)(idx))[3]; \
		idx6 = index; \
		idx7 = index>>32; \
		x1 = (__m128)_mm_cvtsi32_si128(*(int *)(&base[idx0])); \
		x1 = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(x1), *((int *)&base[idx1]),1)); \
		x1 = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(x1), *((int *)&base[idx2]),2)); \
		x1 = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(x1), *((int *)&base[idx3]),3)); \
		x2 = (__m128)_mm_cvtsi32_si128(*(int *)(&base[idx4])); \
		x2 = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(x2), *((int *)&base[idx5]),1)); \
		x2 = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(x2), *((int *)&base[idx6]),2)); \
		x2 = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(x2), *((int *)&base[idx7]),3)); \
		output = _mm256_insertf128_ps(_mm256_castps128_ps256(x1),x2,1); \
	}

// Added by SJP: May not work...
// Certainly won't be "fast" as claimed.
#define fast_scatter_ps(base,idx,data) \
	{ \
		__m128 x1, x2; \
		unsigned long long index; \
		unsigned idx0,idx2,idx4,idx6; \
		unsigned long long idx1,idx3,idx5,idx7; \
	\
		index = ((unsigned long long *)(idx))[0]; \
		idx0 = index; \
		idx1 = index>>32; \
		index = ((unsigned long long *)(idx))[1]; \
		idx2 = index; \
		idx3 = index>>32; \
		index = ((unsigned long long *)(idx))[2]; \
		idx4 = index; \
		idx5 = index>>32; \
		index = ((unsigned long long *)(idx))[3]; \
		idx6 = index; \
		idx7 = index>>32; \
		x1 = _mm256_castps256_ps128(data); \
		base[idx0] = _mm_cvtss_f32(x1); \
		*((int *)&base[idx1]) = _mm_extract_epi32((__m128i) x1, 1); \
		*((int *)&base[idx2]) = _mm_extract_epi32((__m128i)x1, 2); \
		*((int *)&base[idx3]) = _mm_extract_epi32((__m128i)x1, 3); \
		x2 = _mm256_extractf128_ps(data, 1); \
		base[idx4] = _mm_cvtss_f32(x2); \
		*((int *)&base[idx5]) = _mm_extract_epi32((__m128i)x2, 1); \
		*((int *)&base[idx6]) = _mm_extract_epi32((__m128i)x2, 2); \
		*((int *)&base[idx7]) = _mm_extract_epi32((__m128i)x2, 3); \
	}

  fptype inline _MM_REDUCE (SIMDFPTYPE val)
  {
  	  float retval;
	  __m128 xmm0, xmm1;
      xmm0 = _mm256_extractf128_ps(val, 1);
      xmm1 = _mm_add_ps(_mm256_castps256_ps128(val), xmm0);
      xmm0 = _mm_movehl_ps(xmm0, xmm1);
      xmm1 = _mm_add_ps(xmm0, xmm1);
      xmm0 = _mm_movehdup_ps(xmm1);
      xmm1 = _mm_add_ps(xmm0, xmm1);
      _mm_store_ss(&retval, xmm1);
	  return retval;
  }

#elif defined(PRECISION) && (PRECISION == 2)

  #define _masktype unsigned __int64 
  #define fptype double
  #define fpconst(n) (n)  
  #define MPI_FPTYPE MPI_DOUBLE  
  #define VLEN 4
  //#define NNZ_BITS (double)(8*8) 
  #define SIMDFPTYPE __m256d
 
  #undef CEIL_SIMD
  #define CEIL_SIMD(A)  ((int) ceil((double)(A)/(double)VLEN)*VLEN)
 
  #define _MM_SETZERO _mm256_setzero_pd
  #define _MM_SET1 _mm256_set1_pd
  #define _MM_BLENDV _mm256_blendv_pd
  #define _MM_BLEND _mm256_blend_pd
  #define _MM_LOAD _mm256_load_pd
  #define _MM_LOAD1 _mm256_broadcast_sd
  #define _MM_LOADU _mm256_loadu_pd
  #define _MM_STORE _mm256_store_pd
  #define _MM_STORE1 _mm256_store_sd
  #define _MM_STOREU _mm256_storeu_pd
  #define _MM_LOADU _mm256_loadu_pd
 
  //#define _MM_RSQRT _mm256_rsqrt_pd
  #define _MM_RSQRT(src) _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_sqrt_pd(src))
  #define _MM_SQRT _mm256_sqrt_pd
  #define _MM_FLOOR _mm256_floor_pd
  #define _MM_MUL _mm256_mul_pd
  #define _MM_ADD _mm256_add_pd
  #define _MM_SUB _mm256_sub_pd
  #define _MM_MAX _mm256_max_pd
  #define _MM_MIN _mm256_min_pd
  #define _MM_AND _mm256_and_pd
  #define _MM_ANDNOT _mm256_andnot_pd
  #define _MM_OR _mm256_or_pd
  #define _MM_XOR _mm256_xor_pd
  #define _MM_RCP _mm256_rcp_pd
  #define _MM_DIV _mm256_div_pd
  #define _MM_EXP _mm256_exp_pd
  #define _MM_LOG _mm256_log_pd
  #define _MM_CMPEQ(a, b) _mm256_cmp_pd(a, b, _CMP_EQ_OQ)
  #define _MM_CMPNEQ(a, b) _mm256_cmp_pd(a, b, _CMP_NEQ_OQ)
  #define _MM_CMPLT(a, b) _mm256_cmp_pd(a, b, _CMP_LT_OS)
  #define _MM_CMPNLE(a, b) _mm256_cmp_pd(a, b, _CMP_NLE_US)
  #define _MM_CMPNLT(a, b) _mm256_cmp_pd(a, b, _CMP_NLT_US)
  #define _MM_HADD _mm256_hadd_pd
  #define _MM_DP _mm256_dp_pd
  #define _MM_BROADCAST _mm256_broadcast_sd
//  extern "C" { __m256 vmldExp2(__m256 a); }

#ifdef AVX2
  #define _MM_MADD213 	_mm256_fmadd_pd
  #define _MM_MSUB213 	_mm256_fmsub_pd
  #define _MM_MADD132(a, dest, b) 	_mm256_fmadd_pd(a, b, dest)		
  #define _MM_MSUB132(a, dest, b) 	_mm256_fmsub_pd(a, b, dest)		
  #define _MM_MADD231(dest, a, b) 	_mm256_fmadd_pd(a, b, dest)
  #define _MM_MSUB231(dest, a, b) 	_mm256_fmsub_pd(a, b, dest)	// a*b-c
#else
  #define _MM_MADD213(a,b,c)	_MM_ADD(c, _MM_MUL(a, b))
  #define _MM_MSUB213(a,b,c)	_MM_SUB(_MM_MUL(a, b),c)
  #define _MM_MADD132(a, b, c) 	_MM_ADD(b, _MM_MUL(a, c))
  #define _MM_MSUB132(a, b, c) 	_MM_SUB(_MM_MUL(a, c), b)
  #define _MM_MADD231(a, b, c) 	_MM_ADD(a, _MM_MUL(b, c))
  #define _MM_MSUB231(a, b, c) 	_MM_SUB(_MM_MUL(b, c), a)
  #define _MM_MSUBR231(a, b, c)	_MM_SUB(a, _MM_MUL(b, c))
#endif

   #define _MM_MOVEMASK _mm256_movemask_pd
  #define _MM_PERMUTE2 _mm256_permute2_pd
  #define _MM_PERMUTE _mm256_permute_pd
  #define _MM_PERMUTE2F128 _mm256_permute2f128_pd
  #define _MM_INSERTF128 _mm256_insertf128_pd
  #define _MM_EXTRACTF128 _mm256_extractf128_pd
  #define _MM_SHUFPD _mm256_shuffle_pd
  #define _MM_CASTSI_PD _mm256_castsi256_pd
  #define _MM_CASTPD_SI _mm256_castpd_si256
  #define _MM_CAST128TO256 _mm256_castpd128_pd256
  #define _MM_CAST256TO128 _mm256_castpd256_pd128
  #define _MM_UNPACK_LO _mm256_unpacklo_pd
  #define _MM_UNPACK_HI _mm256_unpackhi_pd

  #define ABS fabs
  #define _MM_GATHER _mm256_gather_pd_nomask
  #define _MM_GATHER_Scalar(base, ind, dense, vdense) dense[0]=base[ind[0]]; dense[1]=base[ind[1]]; dense[2]=base[ind[2]]; dense[3]=base[ind[3]];\
                                                      vdense = _MM_LOAD(dense);

  #define MONE -1.0
  #define ZERO 0.0
  #define ONE 1.0
  #define TWO 2.0
  #define THRESH 0.000001

  __declspec(align(64)) static unsigned __int64  MASK[5][4]={{0,0,0,0}, {0x8000000000000000,0,0,0}, 
                                            {0x8000000000000000,0x8000000000000000,0,0}, 
                                            {0x8000000000000000,0x8000000000000000,0x8000000000000000,0}, 
                                            {0x8000000000000000,0x8000000000000000,0x8000000000000000,0x8000000000000000}}; 

	#define fast_gather_pd(base,idx,output) \
	{ \
		__m128 x1, x2; \
		unsigned long long index; \
		unsigned idx0,idx2; \
		unsigned long long idx1,idx3; \
	\
		index = ((unsigned long long *)(idx))[0]; \
		idx0 = index; \
		idx1 = index>>32; \
		index = ((unsigned long long *)(idx))[1]; \
		idx2 = index; \
		idx3 = index>>32; \
		x1 = (__m128)_mm_cvtsi64_si128(*(__int64 *)(&base[idx0])); \
		x1 = _mm_castsi128_ps(_mm_insert_epi64(_mm_castps_si128(x1), *((int *)&base[idx1]),1)); \
		x2 = (__m128)_mm_cvtsi64_si128(*(__int64 *)&base[idx2]); \
		x2 = _mm_castsi128_ps(_mm_insert_epi64(_mm_castps_si128(x2), *((int *)&base[idx3]),1)); \
		output = _mm256_castps_pd(_mm256_insertf128_ps(_mm256_castps128_ps256(x1),x2,1)); \
	}

#else
#  error "Undefined PRECISION: should be 1 (single)  or 2 (double)" 
#endif

#define _MM_FMA _MM_MADD231
#define _MM_FMSUB _MM_MSUB231

#ifndef AVX2
//#define SIMDTYPE_I32 SIMDFPTYPE
#endif
#define SIMDMASK SIMDFPTYPE
//---------------------
// Common for all
//#define MASKTYPE SIMDINTTYPE
//#define _MM_MASKNOT	_MM_ANDNOT
//#define _MM_MASK	_MM_AND
//#define _MM_MERGE	_MM_OR
#define _MM_MASK_MUL (a, mask, b) _MM_MERGE(_MM_MASKNOT(mask, a), _MM_MASK(mask, _MM_MUL (a, b)))
#define _MM_MASK_SUB (a, mask, b) _MM_MERGE(_MM_MASKNOT(mask, a), _MM_MASK(mask, _MM_SUB (a, b)))

#define CSR_COLIND_BITS (double)(8*4)

//#include <gather_impl.h>

/*
#define LABEL(val) label ## val
#define sLABEL(val) LABEL(val)
#define atomic_add(src1, src2) \
    { \ 
    fptype   *accumulator=src1,val=src2; \
__asm                mov             rdi,[accumulator] \
__asm                movd            xmm1, val \
__asm                mov             eax,[rdi]   \
__asm sLABEL(__LINE__): movd         xmm0, eax \
__asm                addss          xmm0, xmm1 \
__asm                movd            ebx, xmm0 \
__asm                lock cmpxchg    [rdi], ebx \
__asm                jnz             sLABEL(__LINE__)\
    } 
#define LABEL(val) label ## val 
#define sLABEL(val) LABEL(val) 
#define atomic_add(src1, src2) \ 
    { \  
    fptype   *accumulator=src1,val=src2; \
__asm                mov             rdi,[accumulator] \
__asm                movd            xmm1, val \
__asm                mov             rax,[rdi]   \
__asm sLABEL(__LINE__): movd         xmm0, rax \
__asm                addsd          xmm0, xmm1 \ 
__asm                movd            rbx, xmm0 \
__asm                lock cmpxchg    [rdi], rbx \
__asm                jnz             sLABEL(__LINE__)\
    }

*/




static void printv_ps(__m256 v, char *str)
{
  int i;
  __declspec(align(64)) float tmp[8];
  printf("%s:", str);
  _mm256_store_ps(tmp, v);
  for(i=0; i < 8; i++)
    printf("[%d]=%f ", i, tmp[i]);
  printf("\n");
 
}
static void printv_pd(__m256d v, char *str)
{
  int i;
  __declspec(align(64)) double tmp[4];
  printf("%s:", str);
  _mm256_store_pd(tmp, v);
  for(i=0; i < 4; i++)
    printf("[%d]=%lf ", i, tmp[i]);
  printf("\n");
}

static void printv_pi(__m256i v, char *str)
{
  int i;
  __declspec(align(64)) int tmp[8];
  printf("%s:", str);
  _mm256_store_si256((__m256i*)tmp, v);
  for(i=0; i < 8; i++)
    printf("[%d]=%x ", i, tmp[i]);
  printf("\n");
}





#endif

