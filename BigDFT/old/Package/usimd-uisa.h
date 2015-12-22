#ifndef __USIMD_512
#define __USIMD_512

#include <assert.h>

//#pragma message "Including usumd-uisa.h"

#ifdef __cplusplus
//extern "C" {
#endif

#include <stdio.h>
#ifdef LRB1
  #ifndef __USE_OLD_NAMES__
    #define __USING_OLD_NAMES__ 0
  #endif
  #include <lmmintrin.h>
  #include "lmmintrin_math.h"
#else
  #include <immintrin.h>
#endif

#if defined LRB2 || defined LRB1
  #define _MM_PREFETCH1(A) _mm_vprefetch1((char *)(A), _MM_PFHINT_NONE);
  #define _MM_PREFETCH2(A) _mm_vprefetch2((char *)(A), _MM_PFHINT_NONE);
#else
  //#define _mm_vprefetch1(A, B)
  //#define _mm_vprefetch2(A, B)
  #define _MM_PREFETCH1(A) _mm_prefetch((char *)(A), _MM_HINT_T0);
  #define _MM_PREFETCH2(A) _mm_prefetch((char *)(A), _MM_HINT_T1);
//  #define _MM_PREFETCH1(A) 
//  #define _MM_PREFETCH2(A) 
#endif

#define MASKTYPE __mmask

#if defined LRB1 || defined LRB2
//  #define _MM_LOAD(ptr) _mm512_loadd((char *)(ptr),_MM_FULLUPC_NONE,_MM_BROADCAST32_NONE,_MM_HINT_NONE)
  #define _MM_SPLAT0(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_AAAA)
  #define _MM_SPLAT1(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_BBBB)
  #define _MM_SPLAT2(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_CCCC)
  #define _MM_SPLAT3(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_DDDD)
  #define _MM_LOADU(reg, mt)\
        reg=_mm512_loadunpackld(reg, (int*)(mt), _MM_FULLUPC_NONE, _MM_HINT_NONE);\
        reg=_mm512_loadunpackhd(reg, (int*)(mt)+16, _MM_FULLUPC_NONE, _MM_HINT_NONE);
  #define _MM_MASK_LOADU(reg, mask, mt) \
		reg=_mm512_mask_loadunpackld(reg, mask, (int*)(mt), _MM_FULLUPC_NONE, _MM_HINT_NONE);\
		reg=_mm512_mask_loadunpackhd(reg, mask, (int*)(mt)+16, _MM_FULLUPC_NONE, _MM_HINT_NONE);
  #define _MM_LOADU_L(DEST, SRC)	DEST = _mm512_loadunpackld (DEST, (int *) SRC, _MM_FULLUPC_NONE, _MM_HINT_NONE);

  #define _MM_STOREU(ptr, reg) \
    	_mm512_packstoreld(ptr,			   reg,_MM_DOWNC_NONE,_MM_HINT_NONE); \
	    _mm512_packstorehd((int *)(ptr)+16,reg,_MM_DOWNC_NONE,_MM_HINT_NONE);
  #define _MM_MASK_STOREU(ptr, mask, reg) \
    	_mm512_mask_packstoreld(ptr,		    mask,reg,_MM_DOWNC_NONE,_MM_HINT_NONE); \
	    _mm512_mask_packstorehd((int *)(ptr)+16,mask,reg,_MM_DOWNC_NONE,_MM_HINT_NONE);

  #ifdef LRB1
      #define _MM_STORE(ptr,reg) _mm512_stored((__m512*)(ptr), (reg),_MM_DOWNC_NONE,_MM_SUBSET32_16,_MM_HINT_NONE)
      #define _MM_MASK_STORE(ptr,mask, reg) _mm512_mask_stored((__m512*)(ptr), mask, (reg),_MM_DOWNC_NONE,_MM_SUBSET32_16,_MM_HINT_NONE)
	  //#define __m512i __m512
	  //#define __m512d __m512
	  #define _mm512_store_pd	 _MM_STORE
	  #define _mm512_extstore_ps	 _MM_STORE
	  #define _mm512_mask_gmaxabs_ps _mm512_mask_maxabs_ps
	  #define _MM_UPCONV_PS_NONE _MM_FULLUPC_NONE
	  #define _MM_UPCONV_EPI32_NONE _MM_FULLUPC_NONE
	  #define _MM_UPCONV_PD_NONE _MM_FULLUPC_NONE
      #define _MM_MASK_GATHER_FP _mm512_mask_gatherd
	  #define _MM_CAST_TO_SI
      #define _MM_CAST_FROM_SI 
      #define _MM_CAST_FROM_PS 
      #define _MM_CAST_FROM_PD 
  #else
    #define _MM_STORE_SI(ptr, reg) _mm512_extstore_epi32(ptr, (reg),_MM_DOWNCONV_EPI32_NONE,_MM_HINT_NONE)
    #define _MM_LOAD_SI(A) _mm512_extload_epi32(A, _MM_UPCONV_EPI32_NONE,_MM_BROADCAST32_NONE,_MM_HINT_NONE)
  #endif
  #define _MM_GATHER_SI _mm512_gatherd
  #define _MM_MUL_SI _mm512_mullo_epi32
  #define _MM_ADD_SI _mm512_add_pi
  #define _MM_FMADD_SI _mm512_madd231_pi
  #define _MM_COUNTBITS_16 _mm_countbits_16
  #define _MM_MADD213_SI16 _mm512_fmadd_epi16
#else
  #define _MM_GATHER_SI _mm512_i32gather_epi32

  #define _MM_LOADU_SI(A) _mm512_load_epi32(A, _MM_UPCONV_EPI32_NONE,_MM_BROADCAST32_NONE,_MM_HINT_NONE)
  #define _MM_LOADU_L(DEST, SRC) _mm512_load_epi32(SRC, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE);

  #define _MM_STOREU_SI _MM_STORE_SI
  #define _MM_SETZERO_SI _mm256_setzero_si512
  #define _MM_MUL_SI _mm512_mullo_epi32
  #define _MM_ADD_SI _mm512_add_epi32
  #define _MM_MUL_SI16 _mm512_mullo_epi16
  #define _MM_ADD_SI16 _mm512_add_epi16
  #define _MM_LOAD1_SI(A) _mm256_castps_si512(_mm512_broadcast_ss(A))
  #define _MM_CMPEQ_SI _mm512_cmpeq_epi32
  #define _MM_CMPGT_SI _mm512_cmpgt_epi32
  #define _MM_OR_SI _mm512_or_si256
  #define _MM_FMADD_SI(v1, v2, v3)	_MM_ADD_SI(v1, _MM_MUL_SI(v2,v3))
  #define _MM_MADD213_SI16(v1, v2, v3)	_MM_ADD_SI16(_MM_MUL_SI16(v1, v2), v3)

  #define _MM_MASK_LOAD_SI(reg, mask, mt) _mm512_mask_load_epi32(reg, mask, mt, _MM_UPCONV_EPI32_NONE,_MM_BROADCAST32_NONE, _MM_HINT_NONE)
  #define _MM_MASK_LOADU_SI(reg, mask, addr) reg=_MM_MASK_LOAD_SI(reg, mask, addr)


  #define _MM_MASK_STORE_SI(ptr, mask, reg) _mm512_mask_store_epi32((__m512*)(ptr), mask, (reg), _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE)
  #define _MM_MASK_STOREU_SI _MM_MASK_STORE_SI

  #define _MM_MASK_GATHER_SI _mm512_mask_i32gather_epi32
  inline int _MM_COUNTBITS_16 ( int i) { return _mm_popcnt_u32(i); }

  #define _MM_SPLAT0(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_AAAA)
  #define _MM_SPLAT1(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_BBBB)
  #define _MM_SPLAT2(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_CCCC)
  #define _MM_SPLAT3(reg) _mm512_swizzle_ps (reg, _MM_SWIZ_REG_DDDD)
  #define _MM_LOAD4(addr) _mm512_extload_ps(addr, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)

#endif

#if defined LRB1
  #define _MM_STORE_SI _MM_STORE
  #define _MM_LOAD_SI  _MM_LOAD
  #define _mm512_load_epi32  _mm512_loadd
  #define _mm512_extstore_epi32 _mm512_stored
#else
  #define _MM_LOAD_SI(A) _mm512_extload_epi32(A, _MM_UPCONV_EPI32_NONE,_MM_BROADCAST32_NONE,_MM_HINT_NONE)
  #define _MM_STORE_SI(ptr, reg) _mm512_extstore_epi32((__m512*)(ptr), (reg),_MM_DOWNCONV_EPI32_NONE,_MM_HINT_NONE)
#endif

  #define ZERO 0.0f
  #define ONE 1.0f
  #define SIMDMASK  	__mmask
  #define SIMDTYPE_I32  __m512i

  #define _MM_SET1_SI _mm512_set_1to16_pi

//#define SIMDINTTYPE __m512
//---------------------------------------------------
// Single precision
//---------------------------------------------------
#if defined(PRECISION) && (PRECISION == 1)
  #define fptype float
  #define fpconst(n) (n ## f)  
  #define MPI_FPTYPE MPI_FLOAT  
  #define VLEN 16
//  #define NNZ_BITS (double)(4*8)
  #define SIMDFPTYPE __m512
//  #define SIMDMASK  SIMDFPTYPE

//  #define vRngGaussian vsRngGaussian
//  #define vExp vsExp
//  #define trmm strmm
//  #define vInv vsInv
 
  #define _MM_SETZERO _mm512_setzero_ps
  #define _MM_SET1 _mm512_set_1to16_ps

//  #define _MM_BLENDV _mm512_blendv_ps
//  #define _MM_BLEND _mm512_blend_ps
  #define _MM_LOAD1(A) _MM_SET1(*(A))
  #define _MM_RSQRT _mm512_rsqrt_ps
  #define _MM_SQRT _mm512_sqrt_ps
//  SIMDFPTYPE _MM_ABS(SIMDFPTYPE A) { __declspec(align(64)) float temp[16]; _mm512_store_ps(temp, A, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); for (int i=0;i<16;i++) temp[i]=fabs(temp[i]); return _mm512_load_ps((char *) (temp), _MM_UPCONV_PS_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE);}
  //#define _MM_ABS(A) _mm512_castsi512_ps(_mm512_absdiff_epi32(_mm512_castps_si512(A),_mm512_castps_si512(_mm512_setzero_ps())))
  #define _MM_ABS(A) _mm512_range_ps(A, A, 0x0a)
  #define _MM_ABSDIFF _mm512_absdiff_ps
// #define _MM_ABS(A) _mm512_gmaxabs_ps(A,A)
#if defined LRB1 || defined LRB2
  #define _MM_FLOOR(a) _mm512_cvt_ps2pi( a, _MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE)
  #define _MM_MASK_FLOOR( old, mask, src) \
					_mm512_mask_cvt_ps2pi(old, mask, src, _MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE)
#else
  #define _MM_FLOOR(a) _mm512_cvt_roundps_epi32( a, _MM_ROUND_MODE_NEAREST)
  #define _MM_MASK_FLOOR( old, mask, src) _mm512_mask_cvt_roundps_epi32(old, mask, src, _MM_ROUND_MODE_NEAREST)
#endif
  
#if defined LRB1 || defined LRB2
  #define _MM_CVT_SI_FP(a) _mm512_cvt_pi2ps( a, _MM_EXPADJ_NONE)
#else
  #define _MM_CVT_SI_FP(a) _mm512_cvt_roundepi32_ps( a, _MM_ROUND_MODE_NEAREST)
#endif

  #define _MM_MUL _mm512_mul_ps
  #define _MM_ADD _mm512_add_ps
  #define _MM_SUB _mm512_sub_ps
  #define _MM_MAX _mm512_max_ps
  #define _MM_MIN _mm512_min_ps

//  #define _MM_AND _mm512_and_ps		Not implemented
//  #define _MM_ANDNOT _mm512_andnot_ps
//  #define _MM_OR _mm512_or_ps
//  #define _MM_XOR _mm512_xor_ps
  #define _MM_DIV _mm512_div_ps
#if defined __USE_SMVL__
  #define _MM_EXP __svml_expf16
  #define _m512_exp_ps __svml_expf16
  extern __m512 __svml_expf16 (__m512);
  #define _m512_log_ps __svml_exp8
  #define _MM_LOG __svml_logf16
  extern __m512 __svml_logf16 (__m512);
#else
  #define _MM_RCP _mm512_rcp_ps
  #define _MM_EXP _mm512_exp_ps
  #define _MM_LOG _mm512_log_ps
#endif
  #define _MM_RSQRT14 _mm512_rsqrt14_ps
  #define _MM_MASK_OR _mm512_mask_or_ps

#ifdef LRB1
  #define _MM_LOAD_FP_X(ptr, bcst) _mm512_loadd((char *)(ptr),_MM_FULLUPC_NONE,bcst,_MM_HINT_NONE)
  __m512 inline _MM_BROADCAST(fptype *addr) {return _mm512_loadd ((char *) (addr), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);}
  __m512 inline _MM_BROADCAST_V(fptype val) {return _mm512_loadd ((char *) (&val), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);}
  #define _MM_MASK_MOV(old, k, v) _mm512_mask_movd(old, k, v)
  #define _MM_CMPLT _mm512_cmplt_ps
#else
  #define _MM_MASK_GATHER_FP _mm512_mask_i32gather_ps
  #define _MM_LOAD_FP_X(addr, bcst) _mm512_extload_ps((char *) (addr), _MM_UPCONV_PS_NONE, bcst, _MM_HINT_NONE)
  #define _MM_MASK_LOAD_FP_X(reg, addr, mask, bcst) _mm512_mask_load_ps(reg, (char *) (addr), mask, _MM_UPCONV_PS_NONE, bcst, _MM_HINT_NONE)
//  #define _MM_BROADCAST(addr) _mm512_load_ps ((char *) (addr), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
  __m512 inline _MM_BROADCAST(const fptype *a) { return _mm512_set_1to16_ps(*a); }
  __m512 inline _MM_BROADCAST_V(const float val) { return _mm512_set_1to16_ps(val); }

  #define _MM_STORE(dst_mem, src_reg) _mm512_store_ps(dst_mem, src_reg);
  #define _MM_MASK_STORE(dst_mem, mask, src_reg) _mm512_mask_store_ps(dst_mem, mask, src_reg);
  #define _MM_MASK_MOV(old, k, v) _mm512_mask_mov_ps(old, k, v)

  #define _MM_LOADU(reg, mt) \
        reg=_mm512_loadunpacklo_ps(reg, (int*)(mt));\
        reg=_mm512_loadunpackhi_ps(reg, (int*)(mt)+16);

  #define _MM_CAST_TO_SI _mm512_castps_si512
  #define _MM_CAST_FROM_SI _mm512_castsi512_ps
  #define _MM_CAST_FROM_PS 
  #define _MM_CAST_FROM_PD _mm512_castpd_ps

  #define _mm512_msubr23c1_ps(a,b) _MM_SUB(_MM_SET1(1.0), _MM_MUL(a,b))
  #define _mm512_mask_msubr23c1_ps(v1, k1, v2, v3) _mm512_mask_sub_ps(v1, k1, _MM_SET1(1.0), _MM_MUL(v2,v3) )

  #define _MM_CMPEQ(a, b) _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ)
  #define _MM_CMPLT(a, b) _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
  #define _MM_CMPNLE(a, b) _mm512_cmp_ps_mask(a, b, _CMP_NLE_US)
#endif
  #define _MM_MADD132 	_mm512_madd132_ps
  #define _MM_MSUB132 	_mm512_msub132_ps

#ifdef AVX3
  #define _MM_MADD231(a, b, c) 	_mm512_fmadd_ps(b, c, a)
#else
  #define _MM_MADD231 	_mm512_madd231_ps
#endif

  #define _MM_MADD213   _mm512_madd213_ps
  #define _MM_MSUB213   _mm512_msub213_ps

  #define _MM_LOAD_FP(addr) _MM_LOAD_FP_X( (addr), _MM_BROADCAST32_NONE)
  #define _MM_MASK_LOAD_FP(reg, mask, addr) _MM_MASK_LOAD_FP_X( reg, mask, addr, _MM_BROADCAST32_NONE)

  #define _MM_REDUCE _mm512_reduce_add_ps


//  #define _MM_HADD _mm512_hadd_ps
//  #define _MM_DP _mm512_dp_ps
//  #define _MM_BLENDV _mm512_blendv_ps
//  extern "C" { __m512 vmlsExp4(__m512 a); }
  

/*#define ABS fabsf
  #define _MM_GATHER _mm512_gather_ps_nomask
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

*/
//---------------------------------------------------
// Double precision
//---------------------------------------------------
#elif defined(PRECISION) && (PRECISION == 2)
  #define _masktype unsigned __int64 
  #define fptype double
  #define fpconst(n) (n)  
  #define MPI_FPTYPE MPI_DOUBLE  
  #define VLEN 8
  //#define NNZ_BITS (double)(8*8) 
  #define SIMDFPTYPE __m512d
 
 
  #define _MM_SETZERO _mm512_setzero_pd
  #define _MM_SET1 _mm512_set_1to8_pd
//  #define _MM_SET1 _mm512_set1_pd
//  #define _MM_BLENDV _mm512_blendv_pd
//  #define _MM_BLEND _mm512_blend_pd

//  #define _MM_BLENDV(v1,v2,mask) _mm512_castsi512_pd(_mm512_mask_or_pi( 						\
//		_mm512_castpd_si512(v1), mask, _mm512_castpd_si512(v2), _mm512_castpd_si512(v2) ) )

#ifdef LRB1
  #define _MM_LOAD_FP_X(ptr, bcst) _mm512_loadq((char *)(ptr),_MM_FULLUPC64_NONE,bcst,_MM_HINT_NONE)
  #define _MM_MASK_MOV(old, k, v) _mm512_mask_movq(old, k, v)
  #define _MM_CMPLT _mm512_cmplt_pd
  inline __m512d _MM_BROADCAST(fptype *addr) {return _mm512_loadq ((char *) (addr), _MM_FULLUPC64_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);}
#else
  #define _MM_LOAD_FP_X(addr, bcast) _mm512_extload_pd((char *) (addr), _MM_UPCONV_PD_NONE, bcast, _MM_HINT_NONE)
  #define _MM_MASK_LOAD_FP_X(reg, addr, mask, bcst) _mm512_mask_load_pd(reg, (char *) (addr), mask, _MM_UPCONV_PS_NONE, bcst, _MM_HINT_NONE)
  #define _MM_STORE(dst_mem, src_reg) _mm512_store_pd(dst_mem, src_reg);
  #define _MM_MASK_STORE(dst_mem, mask, src_reg) _mm512_mask_store_pd(dst_mem, mask, src_reg);
  #define _MM_MASK_MOV(old, k, v) _mm512_mask_mov_pd(old, k, v)

  #define _MM_LOADU(reg, mt) \
        reg=_mm512_loadunpacklo_pd(reg, (int*)(mt));\
        reg=_mm512_loadunpackhi_pd(reg, (int*)(mt)+16);

  #define _MM_BROADCAST(addr) _mm512_load_pd ((char *) (addr), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
  #define _mm512_msubr23c1_pd(a,b) _MM_SUB(_MM_SET1(1.0), _MM_MUL(a,b))
  #define _mm512_mask_msubr23c1_pd(v1, k1, v2, v3) _mm512_mask_sub_pd(v1, k1, _MM_SET1(1.0), _MM_MUL(v2,v3) )

  #define _MM_CMPEQ(a, b) _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
  #define _MM_CMPLT(a, b) _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
  #define _MM_CMPNLE(a, b) _mm512_cmp_pd_mask(a, b, _CMP_NLE_US)

#endif
  #define _MM_MADD132 	_mm512_madd132_pd
  #define _MM_MSUB132 	_mm512_msub132_pd

#ifdef AVX3
  #define _MM_MADD231(a, b, c) 	_mm512_fmadd_pd(b, c, a)
#else
  #define _MM_MADD231 	_mm512_madd231_pd
#endif
  #define _MM_MSUB231 	_mm512_msub231_pd

  #define _MM_MADD213   _mm512_madd213_pd
  #define _MM_MSUB213   _mm512_msub213_pd

  #define _MM_LOAD_FP(addr) _MM_LOAD_FP_X( (addr), _MM_BROADCAST64_NONE)
  #define _MM_MASK_LOAD_FP(reg, mask, addr) _MM_MASK_LOAD_FP_X( reg, mask, addr, _MM_BROADCAST64_NONE)


//  #define _MM_STORE1 _mm512_store_sd
//  #define _MM_STOREU _mm512_storeu_pd
//  #define _MM_LOADU _mm512_loadu_pd
 
  //#define _MM_RSQRT _mm512_rsqrt_pd
  #define _MM_RSQRT(src) _mm512_div_pd(SET1(1.0),_mm512_sqrt_pd(src))
  #define _MM_SQRT _mm512_sqrt_pd
//  #define _MM_FLOOR _mm512_floor_pd
  #define _MM_MUL _mm512_mul_pd
  #define _MM_ADD _mm512_add_pd
  #define _MM_SUB _mm512_sub_pd
  #define _MM_MAX _mm512_max_pd

  #define _MM_RCP _mm512_rcp_pd
  #define _MM_DIV _mm512_div_pd

#if !defined(LRB1) && !defined(LRB2)
  #define _MM_CAST_TO_SI _mm512_castpd_si512
  #define _MM_CAST_FROM_SI _mm512_castsi512_pd
  #define _MM_CAST_FROM_PS _mm512_castps_pd
  #define _MM_CAST_FROM_PD 
#endif


#if defined __USE_SMVL__
  #define _MM_EXP __svml_exp8
  #define _m512_exp_pd __svml_exp8
  extern __m512d __svml_exp8 (__m512d);
  #define _m512_log_pd __svml_exp8
  #define _MM_LOG __svml_log8
  extern __m512d __svml_log8 (__m512d);
#else
  #define _MM_RCP _mm512_rcp_pd
  #define _MM_EXP _mm512_exp_pd
  #define _MM_LOG _mm512_log_pd
#endif
  #define _MM_RSQRT14 _mm512_rsqrt14_pd

  #define _MM_HADD _mm512_hadd_pd
  #define _MM_DP _mm512_dp_pd

/*
  #define ABS fabs
  #define _MM_GATHER _mm512_gather_pd_nomask
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

*/
#else
#  error "Undefined PRECISION: should be 1 (single)  or 2 (double)" 
#endif

#define _MM_LOAD _MM_LOAD_FP
#define _MM_FMA _MM_MADD231
#define _MM_FMSUB _MM_MSUB231
#define _MM_LOAD_X _MM_LOAD_FP_X
#define _MM_BLENDV(v1,v2,mask) _MM_MASK_MOV(v1, mask, v2)

//Common stuff
#ifdef LRB1
  #define _MM_MASK_STOREU_SI _MM_MASK_STOREU
  #define _MM_LOAD_SI _MM_LOAD
  #define _MM_MASK_LOAD_SI _MM_MASK_LOAD
  #define _MM_MASK_LOADU_SI _MM_MASK_LOADU
#else
  #define _MM_AND(v1,v2) _MM_CAST_FROM_SI(_mm512_and_epi32(_MM_CAST_TO_SI(v1), _MM_CAST_TO_SI(v2) ) )
  #define _MM_ANDNOT(v1,v2) _MM_CAST_FROM_SI(_mm512_andnot_epi32(_MM_CAST_TO_SI(v1), _MM_CAST_TO_SI(v2) ) )
  #define _MM_OR(v1,v2) _MM_CAST_FROM_SI(_mm512_or_epi32 (_MM_CAST_TO_SI(v1), _MM_CAST_TO_SI(v2) ) )
#ifdef AVX3
  #define _MM_LOADU1arg _MM_LOAD_FP
  #define _MM_ALIGNR_FP(A,B,C) _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(A),_mm512_castps_si512(B),C))
  //#define _MM_LOADU(addr) _MM_LOAD_FP(addr)
  #define _MM_MASK_LOADU(reg, mask, addr) reg=_MM_MASK_LOAD(reg, addr, mask)

  #define _MM_STOREU _MM_STORE
  #define _MM_MASK_STOREU _MM_MASK_STORE
#endif
  #define _MM_LOAD _MM_LOAD_FP
  #define _MM_MASK_LOAD _MM_MASK_LOAD_FP
#endif

#define CSR_COLIND_BITS (double)(8*4)
static void printv_ps(__m512 v, char *str);

static void print_m512(float *p, char *str)
{
#if (PRECISION==1) // otherwise, code needs to change a little
	__m512 m = _mm512_undefined();
//#ifdef LRB1
	_MM_LOADU(m, p);
	printv_ps(m, str);
	_MM_LOADU(m, p + 1);
	printv_ps(m, str);
	_MM_LOADU(m, p + 2);
	printv_ps(m, str);
//#else
	//m = _MM_LOADU(p);
	//printv_ps(m, str);
	//m = _MM_LOADU(p + 1);
	//printv_ps(m, str);
	//m = _MM_LOADU(p + 2);
	//printv_ps(m, str);
//#endif
#endif
}

#ifdef __INTEL_OFFLOAD
__declspec(target(mic))
#endif
static void printv_epi32(__m512i v, char *str)
{
  int i;
  __declspec(align(64)) int tmp[16];
  printf("%s:", str);
  _MM_STORE_SI(tmp, v);
  for(i=0; i < VLEN; i++)
  {
	  tmp[0] = tmp[i];
      printf("[%d]=%d ", i, tmp[0]);
  }
  printf("\n");
}

#ifdef __INTEL_OFFLOAD
__declspec(target(mic))
#endif
static void printv_ps(__m512 v, char *str)
{
  int i;
  __declspec(align(64)) float tmp[16];
  printf("%s:", str);
  _MM_STOREU(tmp, _MM_CAST_FROM_PS(v));
  for(i=0; i < VLEN; i++)
  {
	  tmp[0] = tmp[i];
      printf("[%d]=%f ", i, tmp[0]);
  }
  printf("\n");
}

#ifdef __INTEL_OFFLOAD
__declspec(target(mic))
#endif
static void printv_pd(__m512d v, char *str)
{
  int i;
  __declspec(align(64)) double tmp[8];
  printf("%s:", str);
  _MM_STOREU(tmp, _MM_CAST_FROM_PD(v));
  for(i=0; i < VLEN; i++)
    printf("[%d]=%lf ", i, tmp[i]);
  printf("\n");
}

#ifdef __INTEL_OFFLOAD
__declspec(target(mic))
#endif
static void printv_si(__m512i v, char *str)
{
  int i;
  __declspec(align(64)) int tmp[16];
  printf("%s:", str);
  _MM_STORE(tmp, _MM_CAST_FROM_SI(v));
  for(i=0; i < VLEN; i++)
    printf("[%d]=%x ", i, tmp[i]);
  printf("\n");
}


static float _MM_REDUCE_PS(__m512 mac)
{
	float res;
#ifdef LRB1
	__m512 temp_mac = _mm512_shuf128x32 (mac, _MM_PERM_AADC, _MM_PERM_DCBA);
	mac = _mm512_add_ps (mac, temp_mac);
	temp_mac = _mm512_shuf128x32 (mac, _MM_PERM_AAAB, _MM_PERM_DCBA);
	mac = _mm512_add_ps (mac, temp_mac);
	mac = _mm512_add_ps (mac, _mm512_swizzle_r32(mac, _MM_SWIZ_REG_CDAB));
	mac = _mm512_add_ps (mac, _mm512_swizzle_r32(mac, _MM_SWIZ_REG_BADC));

	_mm512_stored(&res, mac, _MM_DOWNC_NONE, _MM_SUBSET32_1, _MM_HINT_NONE);
#else
  assert(0);
  res = 0.f;
#endif

	return res;
}

#ifdef __cplusplus
//}
#endif

#endif
