#ifndef __USIMD_128
#define __USIMD_128

//#pragma message ("Including usimd-sse.h")
#include <stdio.h>

#ifdef _WIN32
//#include <emmintrin.h>
#include <xmmintrin.h>
#include <intrin.h>
#else
#include <immintrin.h>
#endif

//#define _MM_ABS _mm_abs_epi32
#define _MM_ABS(A) _mm_castsi128_ps(_mm_abs_epi32(_mm_castps_si128(A)))
#define _MM_LOADU_SI _mm_loadu_si128
#define _MM_LOAD_SI(A) _mm_load_si128((__m128i *)A)
#define _MM_STOREU_SI _mm_storeu_si128
#define _MM_STORE_SI(A,B) _mm_store_si128((__m128i *)A,B)
#define _MM_SETZERO_SI _mm_setzero_si128
#define _MM_MUL_SI _mm_mullo_epi32
#define _MM_ADD_SI _mm_add_epi32
#define _MM_LOAD1_SI _mm_load1_epi32
#define _MM_CMPEQ_SI _mm_cmpeq_epi32
#define _MM_CMPGT_SI _mm_cmpgt_epi32
#define _MM_OR_SI _mm_or_si128
#define _MM_CASTINT_FLOAT _mm_castsi128_ps
#define _MM_LOAD1INT(a) _mm_shuffle_epi32((__m128i)_mm_load_ss(a), 0)

#define _mm_vprefetch1(a, b) ;
#define _mm_vprefetch2(a, b) ;

#define SIMD_ALIGN 128

#if defined(PRECISION) && (PRECISION == 1)
  #define VLEN 4
  #define NNZ_BITS (double)(4*8)
  typedef __m128 SIMDFPTYPE;
  #define SIMDMASK  SIMDFPTYPE
  #define fptype float
  #define fpconst(n) (n ## f)  
  #define MPI_FPTYPE MPI_FLOAT  

  #undef CEIL_SIMD
  #define CEIL_SIMD(A)  ((int) ceil((double)(A)/(double)VLEN)*VLEN)

  #define vRngGaussian vsRngGaussian
  #define vExp vsExp
  #define trmm strmm
  #define vInv vsInv
 
  #define _MM_CASTPS_SI _mm_castps_si128
  #define _MM_CASTSI_PS _mm_castsi128_ps

  #define _MM_SETZERO _mm_setzero_ps
  #define _MM_SET1 _mm_set_ps1
  #define _MM_BLENDV _mm_blendv_ps
  #define _MM_BLEND _mm_blend_ps
  #define _MM_LOAD _mm_load_ps
  #define _MM_LOAD1 _mm_load1_ps
  #define _MM_LOAD4 _mm_load_ps
  #define _MM_LOADU1arg _MM_LOADU
  #define _MM_LOADU _mm_loadu_ps
  #define _MM_LOADU_L _MM_LOADU
  #define _MM_LOADH_PI _mm_loadh_pi
  #define _MM_LOADL_PI _mm_loadl_pi
  #define _MM_LOAD_SCALAR _mm_load_ss
  #define _MM_STORE _mm_store_ps
  #define _MM_STORE1 _mm_store_ss
  #define _MM_STOREU _mm_storeu_ps
 
  #define _MM_RSQRT _mm_rsqrt_ps
  #define _MM_SQRT _mm_sqrt_ps
  #define _MM_FLOOR _mm_floor_ps
  #define _MM_MUL _mm_mul_ps
  #define _MM_ADD _mm_add_ps
  #define _MM_SUB _mm_sub_ps
  #define _MM_MAX _mm_max_ps
  #define _MM_MIN _mm_min_ps
  #define _MM_AND _mm_and_ps
  #define _MM_ANDNOT _mm_andnot_ps
  #define _MM_OR _mm_or_ps
  #define _MM_XOR _mm_xor_ps
  #define _MM_RCP _mm_rcp_ps
  #define _MM_DIV _mm_div_ps
  #define _MM_EXP _mm_exp_ps
  #define _MM_LOG _mm_log_ps
  #define _MM_CMPEQ _mm_cmpeq_ps
  #define _MM_CMPNEQ _mm_cmpneq_ps
  #define _MM_CMPLT _mm_cmplt_ps
  #define _MM_CMPLE _mm_cmple_ps
  #define _MM_CMPNLE _mm_cmpnle_ps
  #define _MM_CMPNLT _mm_cmpnlt_ps
  #define _MM_HADD _mm_hadd_ps
  #define _MM_DP _mm_dp_ps

  inline SIMDFPTYPE _MM_BROADCAST(float *a)	    { return _mm_castsi128_ps( _mm_shuffle_epi32(_mm_castps_si128(_mm_load_ss(a)), 0)); }
  inline SIMDFPTYPE _MM_BROADCAST_V(float a)    { return _mm_castsi128_ps( _mm_shuffle_epi32(_mm_castps_si128(_mm_load_ss(&a)), 0)); }

  #define _MM_CAST_TO_SI _mm_castps_si128
  #define _MM_CAST_FROM_SI _mm_castsi128_ps
   
  #define _MM_MOVEMASK _mm_movemask_ps
  #define _MM_PERMUTE _mm_permute_ps
  #define _MM_SHUFPS _mm_shuffle_ps
  #define _MM_UNPACK_LO _mm_unpacklo_ps
  #define _MM_UNPACK_HI _mm_unpackhi_ps
  #define ABS fabsf
  #define _MM_GATHER _mm_gather_ps_nomask
  #define _MM_GATHER_Scalar(base, ind, dense, vdense) dense[0]=base[ind[0]]; dense[1]=base[ind[1]];\
                                                      dense[2]=base[ind[2]]; dense[3]=base[ind[3]]; \
                                                      vdense = _MM_LOAD(dense);
  #define _MM_SPLAT0(xmm_reg) _mm_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(0,0,0,0))
  #define _MM_SPLAT1(xmm_reg) _mm_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(1,1,1,1))
  #define _MM_SPLAT2(xmm_reg) _mm_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(2,2,2,2))
  #define _MM_SPLAT3(xmm_reg) _mm_shuffle_ps(xmm_reg, xmm_reg, _MM_SHUFFLE(3,3,3,3))

  #define _MM_ALIGNR_FP(A,B,C) _mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(A),_mm_castps_si128(B),C*4))

// Fast gather macro.
#define fast_gather_ps(base,idx,output) \
{ \
    __m128 x1, x2, x3; \
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
    output = _mm_load_ss(&base[idx0]); \
    x1 = _mm_load_ss(&base[idx1]); \
    x2 = _mm_load_ss(&base[idx2]); \
    x3 = _mm_load_ss(&base[idx3]); \
    output = _mm_unpacklo_ps(output, x2); \
    x2 = _mm_unpacklo_ps(x1, x3); \
    output = _mm_unpacklo_ps(output, x2); \
}

// Fast scatter macro.
#define fast_scatter_ps(base,idx,data) \
{ \
    __m128 x1, x2, x3; \
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
    base[idx0] = _mm_cvtss_f32(data) ; \
    *((int *)&base[idx1]) = _mm_extract_epi32((__m128i)data, 1); \
    *((int *)&base[idx2]) = _mm_extract_epi32((__m128i)data, 2); \
    *((int *)&base[idx3]) = _mm_extract_epi32((__m128i)data, 3); \
}

  float inline _MM_REDUCE (SIMDFPTYPE val)
  {
	  SIMDFPTYPE v0 = _MM_SET1(0), v1;
	  float retval;
	  v0 = _mm_movehl_ps(v0, val);
	  v0 = _mm_add_ps(v0, val);
      v1 = _mm_movehdup_ps(v0);
      v0 = _mm_add_ps(v0, v1);
	  _mm_store_ss(&retval, v0);
	  return retval;
  }

  #define MONE -1.0f
  #define ZERO 0.0f
  #define ONE 1.0f
  #define TWO 2.0f
  #define THRESH 0.0001f

  __declspec(align(64)) static int MASK[5][4]={{0,0,0,0}, {0x80000000,0,0,0}, {0x80000000, 0x80000000,0,0},
                         {0x80000000, 0x80000000,0x80000000,0}, 
                         {0x80000000, 0x80000000,0x80000000,0x80000000}};

 
#elif defined(PRECISION) && (PRECISION == 2)
  #define _masktype unsigned __int64 
  #define VLEN 2
  #define NNZ_BITS (double)(8*8)
  #define SIMDFPTYPE __m128d
  #define fptype double
  #define fpconst(n) (n)  
  #define MPI_FPTYPE MPI_DOUBLE  
 
  #undef CEIL_SIMD
  #define CEIL_SIMD(A)  ((int) ceil((double)(A)/(double)VLEN)*VLEN)
 
  #define _MM_CASTPD_SI _mm_castpd_si128
  #define _MM_CASTSI_PD _mm_castsi128_pd

  #define _MM_SETZERO _mm_setzero_pd
  #define _MM_SET1 _mm_set1_pd
  #define _MM_BLENDV _mm_blendv_pd
  #define _MM_BLEND _mm_blend_pd
  #define _MM_LOAD _mm_load_pd
  #define _MM_LOAD1 _mm_load1_pd
  #define _MM_LOADU _mm_loadu_pd
  #define _MM_LOADH_PI _mm_loadh_pd
  #define _MM_LOADL_PI _mm_loadl_pd
  #define _MM_STORE _mm_store_pd
  #define _MM_STORE1 _mm_store_sd
  #define _MM_STOREU _mm_storeu_pd
  #define _MM_LOADU _mm_loadu_pd
  #define _MM_LOAD_SCALAR _mm_load_sd
  //#define _MM_RSQRT _mm_rsqrt_pd
  #define _MM_RSQRT(src) _mm_div_pd(_mm_set1_pd(1.0),_mm_sqrt_pd(src))
  #define _MM_SQRT _mm_sqrt_pd
  #define _MM_FLOOR _mm_floor_pd
  #define _MM_MUL _mm_mul_pd
  #define _MM_ADD _mm_add_pd
  #define _MM_SUB _mm_sub_pd
  #define _MM_MAX _mm_max_pd
  #define _MM_MIN _mm_min_pd
  #define _MM_AND _mm_and_pd
  #define _MM_ANDNOT _mm_andnot_pd
  #define _MM_OR _mm_or_pd
  #define _MM_XOR _mm_xor_pd
  #define _MM_RCP _mm_rcp_pd
  #define _MM_DIV _mm_div_pd
  #define _MM_EXP _mm_exp_pd
  #define _MM_LOG _mm_log_pd
  #define _MM_CMPEQ _mm_cmpeq_pd
  #define _MM_CMPNEQ _mm_cmpneq_pd
  #define _MM_CMPLT _mm_cmplt_pd
  #define _MM_CMPNLE _mm_cmpnle_pd
  #define _MM_CMPNLT _mm_cmpnlt_pd
  #define _MM_HADD _mm_hadd_pd
  #define _MM_DP _mm_dp_pd
  #define _MM_BROADCAST _mm_loaddup_pd

  #define _MM_FMA(dest, src1, src2) _mm_add_pd(dest, _mm_mul_pd(src1,src2))
  #define _MM_FMSUB(a, b, c) _mm_sub_pd(_mm_mul_pd(a, b), c)
  #define _MM_MOVEMASK _mm_movemask_pd
  #define _MM_PERMUTE _mm_permute_pd
  #define _MM_SHUFPD _mm_shuffle_pd
  #define _MM_UNPACK_LO _mm_unpacklo_pd
  #define _MM_UNPACK_HI _mm_unpackhi_pd
  #define ABS fabs
  #define _MM_GATHER _mm_gather_pd_nomask
  #define _MM_GATHER_Scalar(base, ind, dense, vdense) dense[0]=base[ind[0]]; dense[1]=base[ind[1]];\
                                                      vdense = _MM_LOAD(dense);
  #define MONE -1.0
  #define ZERO 0.0
  #define ONE 1.0
  #define TWO 2.0
  #define THRESH 0.000001

  __declspec(align(64)) static unsigned __int64  MASK[3][4]={{0,0}, {0x8000000000000000,0}, 
                                            {0x8000000000000000,0x8000000000000000}}; 


#else
#  error "Undefined PRECISION: should be 1 (single)  or 2 (double)" 
#endif

#define SIMDTYPE_I32 __m128i
#define SIMDMASK SIMDFPTYPE

  #define _MM_MADD213(a,b,c)	_MM_ADD(c, _MM_MUL(a, b))
  #define _MM_MSUB213(a,b,c)	_MM_SUB(_MM_MUL(a, b),c)

  #define _MM_MADD132(a, b, c) 	_MM_ADD(b, _MM_MUL(a, c))
  #define _MM_MSUB132(a, b, c)	_MM_SUB(_MM_MUL(a, c), b)

  #define _MM_MSUB231(a, b, c) 	_MM_SUB(_MM_MUL(b, c), a)
  #define _MM_MADD231(a, b, c) 	_MM_ADD(a, _MM_MUL(b, c))
  #define _MM_MSUBR231(a, b, c)	_MM_SUB(a, _MM_MUL(b, c))
  #define _MM_FMA _MM_MADD231
  #define _MM_FMSUB _MM_MSUB231

#define CSR_COLIND_BITS (double)(8*4)

static inline __m128 _mm_gather_ps_nomask(float* base, __m128i regidx)
{
 
    __m128 voutput;
    int idx0,idx1,idx2,idx3;
 
    idx0=_mm_extract_epi32(regidx,0);
    idx1=_mm_extract_epi32(regidx,1);
    idx2=_mm_extract_epi32(regidx,2);
    idx3=_mm_extract_epi32(regidx,3);
 
    voutput = _mm_castsi128_ps(_mm_cvtsi32_si128(*(int *)(&base[idx0])));
    voutput = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(voutput), *((int *)&base[idx1]),1));
    voutput = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(voutput), *((int *)&base[idx2]),2));
    voutput = _mm_castsi128_ps(_mm_insert_epi32(_mm_castps_si128(voutput), *((int *)&base[idx3]),3));
 
    return voutput;
}
/*
static __inline __m128d _mm_gather_pd_nomask(double* base, __m128i regidx)
{
 
    __m128d voutput;
    int idx0,idx1,idx2,idx3;
 
    idx0=_mm_extract_epi32(regidx,0);
    idx1=_mm_extract_epi32(regidx,1);
    voutput = (__m128d)_mm_cvtsi64_si128(*(__int64 *)(&base[idx0]));
    voutput = _mm_castsi128_pd(_mm_insert_epi64(_mm_castpd_si128(voutput), *((__int64 *)&base[idx1]),1));
    return voutput;
}

*/
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




static void printv_ps(__m128 v, char *str)
{
  int i;
  __declspec(align(64)) float tmp[4];
  printf("%s:", str);
  _mm_store_ps(tmp, v);
  for(i=0; i < 4; i++)
    printf("[%d]=%f ", i, tmp[i]);
  printf("\n");
 
}
static void printv_pd(__m128d v, char *str)
{
  int i;
  __declspec(align(64)) double tmp[2];
  printf("%s:", str);
  _mm_store_pd(tmp, v);
  for(i=0; i < 2; i++)
    printf("[%d]=%lf ", i, tmp[i]);
  printf("\n");
}

static void printv_pi(__m128i v, char *str)
{
  int i;
  __declspec(align(64)) int tmp[4];
  printf("%s:", str);
  _mm_store_si128((__m128i *)tmp, v);
  for(i=0; i < 4; i++)
    printf("[%d]=%x ", i, tmp[i]);
  printf("\n");
}

static float _MM_REDUCE_PS(__m128 mac)
{
	float res;

	__m128 xmm0 = _MM_SET1(0);
  xmm0 = _mm_movehl_ps(xmm0, mac);
	mac = _mm_add_ps(xmm0, mac);
	xmm0 = _mm_movehdup_ps(mac);
	mac = _mm_add_ps(xmm0, mac);
	_mm_store_ss(&res, mac);
	return res;
}


#endif

