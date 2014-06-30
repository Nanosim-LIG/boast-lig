/* Fast 10x10 real matrix/complex vector product, using SSE3 instruction set. */
/* Used by elphon_e (module elphon_sp3d5t_2c) when the preprocessor symbol SSE is defined. */
/* NOTES */
/* Last review 10/01/2012 by YMN. */

#ifdef SSE

#include <stdio.h>
#include <emmintrin.h>
#include <pmmintrin.h>

#define Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(j, lda) \
  B = _mm_load_pd(&b[2*j]); \
  A = _mm_loaddup_pd(a+lda*j); \
  C0 = _mm_add_pd(C0, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+1); \
  C1 = _mm_add_pd(C1, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+2); \
  C2 = _mm_add_pd(C2, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+3); \
  C3 = _mm_add_pd(C3, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+4); \
  C4 = _mm_add_pd(C4, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+5); \
  C5 = _mm_add_pd(C5, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+6); \
  C6 = _mm_add_pd(C6, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+7); \
  C7 = _mm_add_pd(C7, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+8); \
  C8 = _mm_add_pd(C8, _mm_mul_pd(A, B)); \
  A = _mm_loaddup_pd(a+lda*j+9); \
  C9 = _mm_add_pd(C9, _mm_mul_pd(A, B))

#define Cj_equal_Sum_i_in_0_9_Aij_Bi(j, lda) \
  A = _mm_loaddup_pd(a+lda*j); \
  C = _mm_mul_pd(A, B0); \
  A = _mm_loaddup_pd(a+lda*j+1); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B1)); \
  A = _mm_loaddup_pd(a+lda*j+2); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B2)); \
  A = _mm_loaddup_pd(a+lda*j+3); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B3)); \
  A = _mm_loaddup_pd(a+lda*j+4); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B4)); \
  A = _mm_loaddup_pd(a+lda*j+5); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B5)); \
  A = _mm_loaddup_pd(a+lda*j+6); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B6)); \
  A = _mm_loaddup_pd(a+lda*j+7); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B7)); \
  A = _mm_loaddup_pd(a+lda*j+8); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B8)); \
  A = _mm_loaddup_pd(a+lda*j+9); \
  C = _mm_add_pd(C, _mm_mul_pd(A, B9)); \
  _mm_store_pd(&c[2*j], C)

/* Fast 10x1 real matrix/complex vector product, using SSE3 instruction set. */

inline void dmatzvec_n_10x1_sse_(const double * const a, const double * const b, double * const c)
{
  __m128d A, B;
  __m128d C0, C1, C2, C3, C4, C5, C6, C7, C8, C9;

  B = _mm_load_pd(&b[0]);
  A = _mm_loaddup_pd(a);
  C0 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+1);
  C1 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+2);
  C2 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+3);
  C3 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+4);
  C4 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+5);
  C5 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+6);
  C6 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+7);
  C7 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+8);
  C8 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+9);
  C9 = _mm_mul_pd(A, B);

  _mm_store_pd(&c[0], C0);
  _mm_store_pd(&c[2], C1);
  _mm_store_pd(&c[4], C2);
  _mm_store_pd(&c[6], C3);
  _mm_store_pd(&c[8], C4);
  _mm_store_pd(&c[10], C5);
  _mm_store_pd(&c[12], C6);
  _mm_store_pd(&c[14], C7);
  _mm_store_pd(&c[16], C8);
  _mm_store_pd(&c[18], C9);
}

/* Fast 10x10 real matrix/complex vector product, using SSE3 instruction set. */

inline void dmatzvec_n_10x10_sse_(const double * const a, const double * const b, double * const c)
{
  __m128d A, B;
  __m128d C0, C1, C2, C3, C4, C5, C6, C7, C8, C9;

  B = _mm_load_pd(&b[0]);
  A = _mm_loaddup_pd(a);
  C0 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+1);
  C1 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+2);
  C2 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+3);
  C3 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+4);
  C4 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+5);
  C5 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+6);
  C6 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+7);
  C7 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+8);
  C8 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+9);
  C9 = _mm_mul_pd(A, B);

  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(1, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(2, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(3, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(4, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(5, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(6, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(7, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(8, 10);
  Ci_equal_Ci_plus_Aij_Bj_for_i_in_0_9(9, 10);

  _mm_store_pd(&c[0], C0);
  _mm_store_pd(&c[2], C1);
  _mm_store_pd(&c[4], C2);
  _mm_store_pd(&c[6], C3);
  _mm_store_pd(&c[8], C4);
  _mm_store_pd(&c[10], C5);
  _mm_store_pd(&c[12], C6);
  _mm_store_pd(&c[14], C7);
  _mm_store_pd(&c[16], C8);
  _mm_store_pd(&c[18], C9);
}

/* Fast 1x10 real matrix/complex vector product, using SSE3 instruction set. */

inline void dmatzvec_n_1x10_sse_(const double * const a, const double * const b, double * const c)
{
  __m128d A, C;
  __m128d B0, B1, B2, B3, B4, B5, B6, B7, B8, B9;

  B0 = _mm_load_pd(&b[0]);
  B1 = _mm_load_pd(&b[2]);
  B2 = _mm_load_pd(&b[4]);
  B3 = _mm_load_pd(&b[6]);
  B4 = _mm_load_pd(&b[8]);
  B5 = _mm_load_pd(&b[10]);
  B6 = _mm_load_pd(&b[12]);
  B7 = _mm_load_pd(&b[14]);
  B8 = _mm_load_pd(&b[16]);
  B9 = _mm_load_pd(&b[18]);

  A = _mm_loaddup_pd(a);
  C = _mm_mul_pd(A, B0);
  A = _mm_loaddup_pd(a+1);
  C = _mm_add_pd(C, _mm_mul_pd(A, B1));
  A = _mm_loaddup_pd(a+2);
  C = _mm_add_pd(C, _mm_mul_pd(A, B2));
  A = _mm_loaddup_pd(a+3);
  C = _mm_add_pd(C, _mm_mul_pd(A, B3));
  A = _mm_loaddup_pd(a+4);
  C = _mm_add_pd(C, _mm_mul_pd(A, B4));
  A = _mm_loaddup_pd(a+5);
  C = _mm_add_pd(C, _mm_mul_pd(A, B5));
  A = _mm_loaddup_pd(a+6);
  C = _mm_add_pd(C, _mm_mul_pd(A, B6));
  A = _mm_loaddup_pd(a+7);
  C = _mm_add_pd(C, _mm_mul_pd(A, B7));
  A = _mm_loaddup_pd(a+8);
  C = _mm_add_pd(C, _mm_mul_pd(A, B8));
  A = _mm_loaddup_pd(a+9);
  C = _mm_add_pd(C, _mm_mul_pd(A, B9));
  _mm_store_pd(&c[0], C);
}

/* Fast 1x10 transpose real matrix/complex vector product, using SSE3 instruction set. */

inline void dmatzvec_t_1x10_sse_(const double * const a, const double * const b, double * const c)
{
  __m128d A, C;
  __m128d B0, B1, B2, B3, B4, B5, B6, B7, B8, B9;

  B0 = _mm_load_pd(&b[0]);
  B1 = _mm_load_pd(&b[2]);
  B2 = _mm_load_pd(&b[4]);
  B3 = _mm_load_pd(&b[6]);
  B4 = _mm_load_pd(&b[8]);
  B5 = _mm_load_pd(&b[10]);
  B6 = _mm_load_pd(&b[12]);
  B7 = _mm_load_pd(&b[14]);
  B8 = _mm_load_pd(&b[16]);
  B9 = _mm_load_pd(&b[18]);

  Cj_equal_Sum_i_in_0_9_Aij_Bi(0, 10);
}

/* Fast 10x10 transpose real matrix/complex vector product, using SSE3 instruction set. */

inline void dmatzvec_t_10x10_sse_(const double * const a, const double * const b, double * const c)
{
  __m128d A, C;
  __m128d B0, B1, B2, B3, B4, B5, B6, B7, B8, B9;

  B0 = _mm_load_pd(&b[0]);
  B1 = _mm_load_pd(&b[2]);
  B2 = _mm_load_pd(&b[4]);
  B3 = _mm_load_pd(&b[6]);
  B4 = _mm_load_pd(&b[8]);
  B5 = _mm_load_pd(&b[10]);
  B6 = _mm_load_pd(&b[12]);
  B7 = _mm_load_pd(&b[14]);
  B8 = _mm_load_pd(&b[16]);
  B9 = _mm_load_pd(&b[18]);

  Cj_equal_Sum_i_in_0_9_Aij_Bi(0, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(1, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(2, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(3, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(4, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(5, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(6, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(7, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(8, 10);
  Cj_equal_Sum_i_in_0_9_Aij_Bi(9, 10);
}

/* Fast 10x1 transpose real matrix/complex vector product, using SSE3 instruction set. */

inline void dmatzvec_t_10x1_sse_(const double * const a, const double * const b, double * const c)
{
  __m128d A, B;
  __m128d C0, C1, C2, C3, C4, C5, C6, C7, C8, C9;

  B = _mm_load_pd(&b[0]);
  A = _mm_loaddup_pd(a);
  C0 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+1);
  C1 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+2);
  C2 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+3);
  C3 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+4);
  C4 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+5);
  C5 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+6);
  C6 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+7);
  C7 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+8);
  C8 = _mm_mul_pd(A, B);
  A = _mm_loaddup_pd(a+9);
  C9 = _mm_mul_pd(A, B);

  _mm_store_pd(&c[0], C0);
  _mm_store_pd(&c[2], C1);
  _mm_store_pd(&c[4], C2);
  _mm_store_pd(&c[6], C3);
  _mm_store_pd(&c[8], C4);
  _mm_store_pd(&c[10], C5);
  _mm_store_pd(&c[12], C6);
  _mm_store_pd(&c[14], C7);
  _mm_store_pd(&c[16], C8);
  _mm_store_pd(&c[18], C9);
}

/* Fast 1x1, 10x10, 1x10 and 10x1 real matrix/complex vector products. */

void dmatzvec_n_10x10_fast_(const int * const n, const int * const m, const double * const a, const double * const b, double * const c)
{
  if (*n == 10)
  {
    if (*m == 10)
    {
      dmatzvec_n_10x10_sse_(a, b, c);
    }
    else if (*m == 1)
    {
      dmatzvec_n_10x1_sse_(a, b, c);
    }
    else
    {
      printf("dmatzvec_n_10x10_fast : Error, invalid (n = %i, m = %i).\n", *n, *m);
      exit(-1);
    }
  }
  else if (*n == 1)
  {
    if (*m == 10)
    {
      dmatzvec_n_1x10_sse_(a, b, c);
    }
    else if (*m == 1)
    {
      c[0] = a[0]*b[0];
      c[1] = a[0]*b[1];
    }
    else
    {
      printf("dmatzvec_n_10x10_fast : Error, invalid (n = %i, m = %i).\n", *n, *m);
      exit(-1);
    }
  }
  else
  {
    printf("dmatzvec_n_10x10_fast : Error, invalid (n = %i, m = %i).\n", *n, *m);
    exit(-1);
  }
}

/* Fast 1x1, 10x10, 1x10 and 10x1 transpose real matrix/complex vector products. */

void dmatzvec_t_10x10_fast_(const int * const n, const int * const m, const double * const a, const double * const b, double * const c)
{
  if (*n == 10)
  {
    if (*m == 10)
    {
      dmatzvec_t_10x10_sse_(a, b, c);
    }
    else if (*m == 1)
    {
      dmatzvec_t_10x1_sse_(a, b, c);
    }
    else
    {
      printf("dmatzvec_t_10x10_fast : Error, invalid (n = %i, m = %i).\n", *n, *m);
      exit(-1);
    }
  }
  else if (*n == 1)
  {
    if (*m == 10)
    {
      dmatzvec_t_1x10_sse_(a, b, c);
    }
    else if (*m == 1)
    {
      c[0] = a[0]*b[0];
      c[1] = a[0]*b[1];
    }
    else
    {
      printf("dmatzvec_t_10x10_fast : Error, invalid (n = %i, m = %i).\n", *n, *m);
      exit(-1);
    }
  }
  else
  {
    printf("dmatzvec_t_10x10_fast : Error, invalid (n = %i, m = %i).\n", *n, *m);
    exit(-1);
  }
}

// /* Fast 10x10 real matrix/complex vector product, using SSE3 instruction set. */
//
// //inline void dmatzvec_n_10x10_sse_(const int * const n, const int * const m, const double * const a, const int * const lda, const double * const b, double * const c)
// inline void dmatzvec_n_10x10_sse_(const double * const a, const double * const b, double * const c)
// {
//   int i;
//   const double * a_ptr;
//   __m128d A, B;
//   __m128d C0, C1, C2, C3, C4, C5, C6, C7, C8, C9;
//
//   a_ptr = a;
//   _mm_prefetch(a_ptr+(10), _MM_HINT_NTA);
//
//   B = _mm_load_pd(&b[0]);
//   A = _mm_loaddup_pd(a_ptr);
//   C0 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+1);
//   C1 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+2);
//   C2 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+3);
//   C3 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+4);
//   C4 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+5);
//   C5 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+6);
//   C6 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+7);
//   C7 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+8);
//   C8 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+9);
//   C9 = _mm_mul_pd(A, B);
//
//   for (i = 1; i < 10; i++)
//   {
//     a_ptr += (10);
//     _mm_prefetch(a_ptr+(10), _MM_HINT_NTA);
//
//     B = _mm_load_pd(&b[2*i]);
//     A = _mm_loaddup_pd(a_ptr);
//     C0 = _mm_add_pd(C0, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+1);
//     C1 = _mm_add_pd(C1, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+2);
//     C2 = _mm_add_pd(C2, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+3);
//     C3 = _mm_add_pd(C3, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+4);
//     C4 = _mm_add_pd(C4, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+5);
//     C5 = _mm_add_pd(C5, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+6);
//     C6 = _mm_add_pd(C6, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+7);
//     C7 = _mm_add_pd(C7, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+8);
//     C8 = _mm_add_pd(C8, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+9);
//     C9 = _mm_add_pd(C9, _mm_mul_pd(A, B));
//   }
//
//   _mm_store_pd(&c[0], C0);
//   _mm_store_pd(&c[2], C1);
//   _mm_store_pd(&c[4], C2);
//   _mm_store_pd(&c[6], C3);
//   _mm_store_pd(&c[8], C4);
//   _mm_store_pd(&c[10], C5);
//   _mm_store_pd(&c[12], C6);
//   _mm_store_pd(&c[14], C7);
//   _mm_store_pd(&c[16], C8);
//   _mm_store_pd(&c[18], C9);
// }
//
// /* Fast 10x10 transpose real matrix/complex vector product, using SSE3 instruction set. */
// /* Columnwise implementation. */
//
// //inline void dmatzvec_t_10x10_sse_1_(const int * const n, const int * const m, const double * const a, const int * const lda, const double * const b, double * const c)
// inline void dmatzvec_t_10x10_sse_1_(const double * const a, const double * const b, double * const c)
// {
//   int i;
//   const double * a_ptr;
//   __m128d A, B;
//   __m128d C0, C1, C2, C3, C4, C5, C6, C7, C8, C9;
//
//   a_ptr = a;
//
//   B = _mm_load_pd(&b[0]);
//   A = _mm_loaddup_pd(a_ptr);
//   C0 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+(10));
//   C1 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+2*(10));
//   C2 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+3*(10));
//   C3 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+4*(10));
//   C4 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+5*(10));
//   C5 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+6*(10));
//   C6 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+7*(10));
//   C7 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+8*(10));
//   C8 = _mm_mul_pd(A, B);
//   A = _mm_loaddup_pd(a_ptr+9*(10));
//   C9 = _mm_mul_pd(A, B);
//
//   for (i = 1; i < 10; i++)
//   {
//     a_ptr++;
//
//     B = _mm_load_pd(&b[2*i]);
//     A = _mm_loaddup_pd(a_ptr);
//     C0 = _mm_add_pd(C0, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+(10));
//     C1 = _mm_add_pd(C1, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+2*(10));
//     C2 = _mm_add_pd(C2, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+3*(10));
//     C3 = _mm_add_pd(C3, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+4*(10));
//     C4 = _mm_add_pd(C4, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+5*(10));
//     C5 = _mm_add_pd(C5, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+6*(10));
//     C6 = _mm_add_pd(C6, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+7*(10));
//     C7 = _mm_add_pd(C7, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+8*(10));
//     C8 = _mm_add_pd(C8, _mm_mul_pd(A, B));
//     A = _mm_loaddup_pd(a_ptr+9*(10));
//     C9 = _mm_add_pd(C9, _mm_mul_pd(A, B));
//   }
//
//   _mm_store_pd(&c[0], C0);
//   _mm_store_pd(&c[2], C1);
//   _mm_store_pd(&c[4], C2);
//   _mm_store_pd(&c[6], C3);
//   _mm_store_pd(&c[8], C4);
//   _mm_store_pd(&c[10], C5);
//   _mm_store_pd(&c[12], C6);
//   _mm_store_pd(&c[14], C7);
//   _mm_store_pd(&c[16], C8);
//   _mm_store_pd(&c[18], C9);
// }
//
// /* Fast 10x10 transpose real matrix/complex vector product, using SSE3 instruction set. */
// /* Linewise implementation. */
//
// //inline void dmatzvec_t_10x10_sse_2_(const int * const n, const int * const m, const double * const a, const int * const lda, const double * const b, double * const c)
// inline void dmatzvec_t_10x10_sse_2_(const double * const a, const double * const b, double * const c)
// {
//   int i;
//   const double * a_ptr;
//   __m128d A, C;
//   __m128d B0, B1, B2, B3, B4, B5, B6, B7, B8, B9;
//
//   a_ptr = a;
//   _mm_prefetch(a_ptr, _MM_HINT_NTA);
//
//   B0 = _mm_load_pd(&b[0]);
//   B1 = _mm_load_pd(&b[2]);
//   B2 = _mm_load_pd(&b[4]);
//   B3 = _mm_load_pd(&b[6]);
//   B4 = _mm_load_pd(&b[8]);
//   B5 = _mm_load_pd(&b[10]);
//   B6 = _mm_load_pd(&b[12]);
//   B7 = _mm_load_pd(&b[14]);
//   B8 = _mm_load_pd(&b[16]);
//   B9 = _mm_load_pd(&b[18]);
//
//   for (i = 0; i < 10; i++)
//   {
//     _mm_prefetch(a_ptr+(10), _MM_HINT_NTA);
//
//     A = _mm_loaddup_pd(a_ptr);
//     C = _mm_mul_pd(A, B0);
//     A = _mm_loaddup_pd(a_ptr+1);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B1));
//     A = _mm_loaddup_pd(a_ptr+2);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B2));
//     A = _mm_loaddup_pd(a_ptr+3);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B3));
//     A = _mm_loaddup_pd(a_ptr+4);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B4));
//     A = _mm_loaddup_pd(a_ptr+5);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B5));
//     A = _mm_loaddup_pd(a_ptr+6);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B6));
//     A = _mm_loaddup_pd(a_ptr+7);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B7));
//     A = _mm_loaddup_pd(a_ptr+8);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B8));
//     A = _mm_loaddup_pd(a_ptr+9);
//     C = _mm_add_pd(C, _mm_mul_pd(A, B9));
//
//     _mm_store_pd(&c[2*i], C);
//
//     a_ptr += (10);
//   }
// }
#endif

