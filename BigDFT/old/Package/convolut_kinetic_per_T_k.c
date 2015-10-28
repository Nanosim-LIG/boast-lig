#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>
#include <boost/preprocessor/repetition/repeat.hpp>
#define __USE_MISC
#include <sys/mman.h>

#ifdef __MIC__
#define PRECISION 2
#include "skbintrin.h"
#endif

#define MULTIPLE_OF_CACHE_LINE(x) ((((x) + 7)/8)*8)
#define IN_STRIDE1(n1) MULTIPLE_OF_CACHE_LINE((n1) + (-lowfil) + (lupfil))
#define IN_STRIDE3(n3) ((n3) + (-lowfil) + (lupfil))
#define IN_STRIDE21(n2, n1) (((n2) + (-lowfil) + (lupfil)) * IN_STRIDE1(n1) + 24)
// Magic offset 24 to avoid conflict misses primarily on KNC
#define OFFSET_OF_INPUT(z, y, x, n2, n1) ((x) * IN_STRIDE21(n2, n1) + (y) * IN_STRIDE1(n1) + z)

#define OUT_STRIDE1(n1) MULTIPLE_OF_CACHE_LINE(n1)
#define OUT_STRIDE2(n2) (n2)
#define OUT_STRIDE3(n3) (n3)
#define OUT_STRIDE21(n2, n1) ((n2) * OUT_STRIDE1(n1))

#define OFFSET_OF_OUTPUT(z, y, x, n2, n1) ((x) * OUT_STRIDE21(n2, n1) + (y) * OUT_STRIDE1(n1) + z)

void barrier(int threadid);

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define L1_PREFETCH_DISTANCE 2

static const double DEFAULT_CPU_FREQ = 3.33e9;
double get_cpu_freq()
{
  static double freq = DBL_MAX;
  if (DBL_MAX == freq) {
    volatile double a = rand()%1024, b = rand()%1024;
    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);
    unsigned long long t1 = __rdtsc();
    for (size_t i = 0; i < 1024L*1024; i++) {
      a += a*b + b/a;
    }
    unsigned long long dt = __rdtsc() - t1;
    gettimeofday(&tv2, NULL);
    freq = dt/((tv2.tv_sec - tv1.tv_sec) + (tv2.tv_usec - tv1.tv_usec)/1.e6);
  }

  return freq;
}

int nthreads;
static unsigned long long time_init = 0;
static unsigned long long time_t1 = 0, time_t2 = 0, time_t3 = 0;

// !>   Applies the modified kinetic energy operator onto x to get y. 
// !!   Computes the kinetic energy too.
// !!   Works for periodic BC.
// !!   Modified kinetic energy operator:
// !!   A=-1/2 exp(Ikr) Delta exp(-Ikr)
// !!   where k=(k1,k2,k3); r=(x,y,z)
// !!
// !! 
// !dir$ attributes offload:mic :: convolut_kinetic_per_T_k
void convolut_kinetic_per_T_k(const unsigned int n1,
                              const unsigned int n2,
                              const unsigned int n3,
                              const double *hgrid,
                              const double *x1,
                              const double *x2,
                              double *y1,
                              double *y2,
                              double *kstrten,
                              const double k1,
                              const double k2,
                              const double k3,
                              unsigned threadid)
{
  unsigned long long time_1, time_2, time_3;
  time_t1 = time_t2 = time_t3 = 0;

  time_1 = __rdtsc();

  // When we integrate these should be defined by the app already.
  const int lowfil=-14;
  const int lupfil=14;

  double scale[3];
  scale[0]=(-.5/(hgrid[0]*hgrid[0]));
  scale[1]=(-.5/(hgrid[1]*hgrid[1]));
  scale[2]=(-.5/(hgrid[2]*hgrid[2]));

  double scale1[3];
  scale1[0]=k1/hgrid[0];
  scale1[1]=k2/hgrid[1];
  scale1[2]=k3/hgrid[2];

  double cx, cy, cz;
  cx=0.5*k1*k1;
  cy=0.5*k2*k2;
  cz=0.5*k3*k3;

  // second derivative filters for Daubechies 16
  const unsigned int offset = -(lowfil);

  double fil1[3][29];
  fil1[0][offset+0]=   -3.5536922899131901941296809374e0*scale[0];
  fil1[1][offset+0]=   -3.5536922899131901941296809374e0*scale[1];
  fil1[2][offset+0]=   -3.5536922899131901941296809374e0*scale[2];

  fil1[0][offset+1]=    2.2191465938911163898794546405e0*scale[0];
  fil1[1][offset+1]=    2.2191465938911163898794546405e0*scale[1];
  fil1[2][offset+1]=    2.2191465938911163898794546405e0*scale[2];

  fil1[0][offset+2]=   -0.6156141465570069496314853949e0*scale[0];
  fil1[1][offset+2]=   -0.6156141465570069496314853949e0*scale[1];
  fil1[2][offset+2]=   -0.6156141465570069496314853949e0*scale[2];

  fil1[0][offset+3]=    0.2371780582153805636239247476e0*scale[0];
  fil1[1][offset+3]=    0.2371780582153805636239247476e0*scale[1];
  fil1[2][offset+3]=    0.2371780582153805636239247476e0*scale[2];

  fil1[0][offset+4]=   -0.0822663999742123340987663521e0*scale[0];
  fil1[1][offset+4]=   -0.0822663999742123340987663521e0*scale[1];
  fil1[2][offset+4]=   -0.0822663999742123340987663521e0*scale[2];

  fil1[0][offset+5]=    0.02207029188482255523789911295638968409e0*scale[0];
  fil1[1][offset+5]=    0.02207029188482255523789911295638968409e0*scale[1];
  fil1[2][offset+5]=    0.02207029188482255523789911295638968409e0*scale[2];

  fil1[0][offset+6]=   -0.409765689342633823899327051188315485e-2*scale[0];
  fil1[1][offset+6]=   -0.409765689342633823899327051188315485e-2*scale[1];
  fil1[2][offset+6]=   -0.409765689342633823899327051188315485e-2*scale[2];

  fil1[0][offset+7]=    0.45167920287502235349480037639758496e-3*scale[0];
  fil1[1][offset+7]=    0.45167920287502235349480037639758496e-3*scale[1];
  fil1[2][offset+7]=    0.45167920287502235349480037639758496e-3*scale[2];

  fil1[0][offset+8]=   -0.2398228524507599670405555359023135e-4*scale[0];
  fil1[1][offset+8]=   -0.2398228524507599670405555359023135e-4*scale[1];
  fil1[2][offset+8]=   -0.2398228524507599670405555359023135e-4*scale[2];

  fil1[0][offset+9]=    2.0904234952920365957922889447361e-6*scale[0];
  fil1[1][offset+9]=    2.0904234952920365957922889447361e-6*scale[1];
  fil1[2][offset+9]=    2.0904234952920365957922889447361e-6*scale[2];

  fil1[0][offset+10]=  -3.7230763047369275848791496973044e-7*scale[0];
  fil1[1][offset+10]=  -3.7230763047369275848791496973044e-7*scale[1];
  fil1[2][offset+10]=  -3.7230763047369275848791496973044e-7*scale[2];

  fil1[0][offset+11]=  -1.05857055496741470373494132287e-8*scale[0];
  fil1[1][offset+11]=  -1.05857055496741470373494132287e-8*scale[1];
  fil1[2][offset+11]=  -1.05857055496741470373494132287e-8*scale[2];

  fil1[0][offset+12]=  -5.813879830282540547959250667e-11*scale[0];
  fil1[1][offset+12]=  -5.813879830282540547959250667e-11*scale[1];
  fil1[2][offset+12]=  -5.813879830282540547959250667e-11*scale[2];

  fil1[0][offset+13]=   2.70800493626319438269856689037647576e-13*scale[0];
  fil1[1][offset+13]=   2.70800493626319438269856689037647576e-13*scale[1];
  fil1[2][offset+13]=   2.70800493626319438269856689037647576e-13*scale[2];

  fil1[0][offset+14]=  -6.924474940639200152025730585882e-18*scale[0];
  fil1[1][offset+14]=  -6.924474940639200152025730585882e-18*scale[1];
  fil1[2][offset+14]=  -6.924474940639200152025730585882e-18*scale[2];

  for (int k = 0; k < (-lowfil); ++k)
  {
      fil1[0][k] = fil1[0][lupfil + (-lowfil) -k];
      fil1[1][k] = fil1[1][lupfil + (-lowfil) -k];
      fil1[2][k] = fil1[2][lupfil + (-lowfil) -k];
  }

  double fil2[3][29];
  fil2[0][offset+0]= 0.;
  fil2[1][offset+0]= 0.;
  fil2[2][offset+0]= 0.;

  fil2[0][offset+1]= 0.8834460460908270942785856e0*scale1[0];
  fil2[1][offset+1]= 0.8834460460908270942785856e0*scale1[1];
  fil2[2][offset+1]= 0.8834460460908270942785856e0*scale1[2];

  fil2[0][offset+2]= -0.3032593514765938346887962e0*scale1[0];
  fil2[1][offset+2]= -0.3032593514765938346887962e0*scale1[1];
  fil2[2][offset+2]= -0.3032593514765938346887962e0*scale1[2];

  fil2[0][offset+3]= 0.1063640682894442760934532e0*scale1[0];
  fil2[1][offset+3]= 0.1063640682894442760934532e0*scale1[1];
  fil2[2][offset+3]= 0.1063640682894442760934532e0*scale1[2];

  fil2[0][offset+4]= -0.03129014783948023634381564e0*scale1[0];
  fil2[1][offset+4]= -0.03129014783948023634381564e0*scale1[1];
  fil2[2][offset+4]= -0.03129014783948023634381564e0*scale1[2];

  fil2[0][offset+5]= 0.006958379116450707495020408e0*scale1[0];
  fil2[1][offset+5]= 0.006958379116450707495020408e0*scale1[1];
  fil2[2][offset+5]= 0.006958379116450707495020408e0*scale1[2];

  fil2[0][offset+6]= -0.001031530213375445369097965e0*scale1[0];
  fil2[1][offset+6]= -0.001031530213375445369097965e0*scale1[1];
  fil2[2][offset+6]= -0.001031530213375445369097965e0*scale1[2];

  fil2[0][offset+7]= 0.00007667706908380351933901775e0*scale1[0];
  fil2[1][offset+7]= 0.00007667706908380351933901775e0*scale1[1];
  fil2[2][offset+7]= 0.00007667706908380351933901775e0*scale1[2];

  fil2[0][offset+8]= 2.451992111053665419191564e-7*scale1[0];
  fil2[1][offset+8]= 2.451992111053665419191564e-7*scale1[1];
  fil2[2][offset+8]= 2.451992111053665419191564e-7*scale1[2];

  fil2[0][offset+9]= 3.993810456408053712133667e-8*scale1[0];
  fil2[1][offset+9]= 3.993810456408053712133667e-8*scale1[1];
  fil2[2][offset+9]= 3.993810456408053712133667e-8*scale1[2];

  fil2[0][offset+10]=-7.207948238588481597101904e-8*scale1[0];
  fil2[1][offset+10]=-7.207948238588481597101904e-8*scale1[1];
  fil2[2][offset+10]=-7.207948238588481597101904e-8*scale1[2];

  fil2[0][offset+11]=-9.697184925637300947553069e-10*scale1[0];
  fil2[1][offset+11]=-9.697184925637300947553069e-10*scale1[1];
  fil2[2][offset+11]=-9.697184925637300947553069e-10*scale1[2];

  fil2[0][offset+12]=-7.252206916665149851135592e-13*scale1[0];
  fil2[1][offset+12]=-7.252206916665149851135592e-13*scale1[1];
  fil2[2][offset+12]=-7.252206916665149851135592e-13*scale1[2];

  fil2[0][offset+13]=1.240078536096648534547439e-14*scale1[0];
  fil2[1][offset+13]=1.240078536096648534547439e-14*scale1[1];
  fil2[2][offset+13]=1.240078536096648534547439e-14*scale1[2];

  fil2[0][offset+14]=-1.585464751677102510097179e-19*scale1[0];
  fil2[1][offset+14]=-1.585464751677102510097179e-19*scale1[1];
  fil2[2][offset+14]=-1.585464751677102510097179e-19*scale1[2];

  for (int k = 0; k < (-lowfil); ++k)
  {
      fil2[0][k] = -fil2[0][lupfil + (-lowfil) -k];
      fil2[1][k] = -fil2[1][lupfil + (-lowfil) -k];
      fil2[2][k] = -fil2[2][lupfil + (-lowfil) -k];
  }

  time_2 = __rdtsc();

  if (0 == threadid)
  time_init += (time_2 - time_1);

  kstrten[0]=0.0;
  kstrten[1]=0.0;
  kstrten[2]=0.0;
  kstrten[3]=0.0;
  kstrten[4]=0.0;
  kstrten[5]=0.0;

  double kstrt1, kstrt2, kstrt3;
  kstrt1 = kstrt2 = kstrt3 = 0.0;

  double tt1, tt2;
  tt1 = tt2 = 0;

  barrier(threadid);
  time_1 = __rdtsc();

  int i2_per_thread = (n2 + nthreads - 1)/nthreads;
  int i2_begin = i2_per_thread*threadid;
  int i2_end = MIN(i2_begin + i2_per_thread, n2);

  for (int i2 = i2_begin; i2 < i2_end; ++i2)
  {
     // (1/2) |d/dx+ik_x)|^2
#define INTRINSIC
#ifdef INTRINSIC
     for (int i1 = 0; i1 < n1; i1 += 8)
#else
     for (int i1 = 0; i1 < n1; i1++)
#endif
     {
         // for (int i1 = 0; i1 < n1; ++i1)
         for (int i3 = 0; i3 < n3; i3+=2)
         // for (int i1 = 0; i1 < n1; i1+=4)
         {
#ifdef INTRINSIC
             size_t base_addrA = OFFSET_OF_INPUT(i1 - lowfil, i2 - lowfil, i3 - lowfil, n2, n1);
             size_t base_addrB = OFFSET_OF_INPUT(i1 - lowfil, i2 - lowfil, i3 + 1 - lowfil, n2, n1);

             SIMDFPTYPE x10A = _MM_LOAD(x1 + base_addrA);
             SIMDFPTYPE x20A = _MM_LOAD(x2 + base_addrA);
             SIMDFPTYPE tt1A = _MM_MUL(x10A, _MM_SET1(cz));
             SIMDFPTYPE tt2A = _MM_MUL(x20A, _MM_SET1(cz));

             SIMDFPTYPE x10B = _MM_LOAD(x1 + base_addrB);
             SIMDFPTYPE x20B = _MM_LOAD(x2 + base_addrB);
             SIMDFPTYPE tt1B = _MM_MUL(x10B, _MM_SET1(cz));
             SIMDFPTYPE tt2B = _MM_MUL(x20B, _MM_SET1(cz));

#define _MM_MSUBR231(a, b, c) _mm512_fnmadd_pd(b, c, a)

             SIMDFPTYPE coeff1, coeff2, x11, x22;
             const size_t in_stride = IN_STRIDE21(n2, n1);

             base_addrA += lowfil*in_stride;
             base_addrB += lowfil*in_stride;

             const size_t out_stride = OUT_STRIDE21(n2, n1);
             size_t o_base_addrA = OFFSET_OF_OUTPUT(i1, i2, i3, n2, n1);
             size_t o_base_addrB = OFFSET_OF_OUTPUT(i1, i2, i3 + 1, n2, n1);

             _MM_PREFETCH1(x1 + base_addrA + (28 + L1_PREFETCH_DISTANCE)*in_stride);
             _MM_PREFETCH1(x2 + base_addrA + (28 + L1_PREFETCH_DISTANCE)*in_stride);
             _MM_PREFETCH1(y1 + o_base_addrA + L1_PREFETCH_DISTANCE*out_stride);
             _MM_PREFETCH1(y2 + o_base_addrA + L1_PREFETCH_DISTANCE*out_stride);


#define MAD_FILTERA(z, l, dim) \
             coeff1 = _MM_SET1(fil1[dim][l]), coeff2 = _MM_SET1(fil2[dim][l]); \
             x11 = _MM_LOAD(x1 + base_addrA + l*in_stride); \
             x22 = _MM_LOAD(x2 + base_addrA + l*in_stride); \
             tt1A = _MM_MSUBR231(_MM_FMA(tt1A, x11, coeff1), x22, coeff2); \
             tt2A = _MM_FMA(_MM_FMA(tt2A, x22, coeff1), x11, coeff2);

#define MAD_FILTERB(z, l, dim) \
             coeff1 = _MM_SET1(fil1[dim][l]), coeff2 = _MM_SET1(fil2[dim][l]); \
             x11 = _MM_LOAD(x1 + base_addrB + l*in_stride); \
             x22 = _MM_LOAD(x2 + base_addrB + l*in_stride); \
             tt1B = _MM_MSUBR231(_MM_FMA(tt1B, x11, coeff1), x22, coeff2); \
             tt2B = _MM_FMA(_MM_FMA(tt2B, x22, coeff1), x11, coeff2);

             BOOST_PP_REPEAT(29, MAD_FILTERA, 2);
             BOOST_PP_REPEAT(29, MAD_FILTERB, 2);

             _mm512_storenrngo_pd(
               y1 + o_base_addrA, _MM_ADD(_MM_LOAD(y1 + o_base_addrA), tt1A));
             _mm512_storenrngo_pd(
               y2 + o_base_addrA, _MM_ADD(_MM_LOAD(y2 + o_base_addrA), tt2A));

             _mm512_storenrngo_pd(
               y1 + o_base_addrB, _MM_ADD(_MM_LOAD(y1 + o_base_addrB), tt1B));
             _mm512_storenrngo_pd(
               y2 + o_base_addrB, _MM_ADD(_MM_LOAD(y2 + o_base_addrB), tt2B));

             SIMDFPTYPE kstrt1_v = _MM_FMA(_MM_MUL(tt1A, x10A), tt2A, x20A);
             kstrt1 += _mm512_reduce_add_pd(kstrt1_v);

             kstrt1_v = _MM_FMA(_MM_MUL(tt1B, x10B), tt2B, x20B);
             kstrt1 += _mm512_reduce_add_pd(kstrt1_v);

#else
// FIXME : This will not work now - have to be fixed
             tt1 = x1[OFFSET_OF_INPUT(i1 + (-lowfil), i2 + (-lowfil), i3 + (-lowfil), n2, n1)] * cx;
             tt2 = x2[OFFSET_OF_INPUT(i1 + (-lowfil), i2 + (-lowfil) ,i3 + (-lowfil), n2, n1)] * cx;
             for (int l = lowfil; l <= lupfil; ++l)
             {
                 int j = (i1 + l);

                 tt1+=x1[OFFSET_OF_INPUT(j + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*fil1[l + (-lowfil)][0]-
                      x2[OFFSET_OF_INPUT(j + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*fil2[l + (-lowfil)][0];
                 tt2+=x2[OFFSET_OF_INPUT(j + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*fil1[l + (-lowfil)][0]+
                      x1[OFFSET_OF_INPUT(j + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*fil2[l + (-lowfil)][0];
             }
             y1[OFFSET_OF_OUTPUT(i1, i2, i3, n2, n1)] += tt1;
             y2[OFFSET_OF_OUTPUT(i1, i2, i3, n2, n1)] += tt2;

             kstrt1+=tt1*x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]+
                     tt2*x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)];
#endif
        }
     }
  }

  barrier(threadid);
  time_2 = __rdtsc();

     // + (1/2) d^2/dy^2
  int i3_per_thread = (n3 + nthreads - 1)/nthreads;
  int i3_begin = i3_per_thread*threadid;
  int i3_end = MIN(i3_begin + i3_per_thread, n3);

  for (int i3 = i3_begin; i3 < i3_end; ++i3)
  {
#ifdef INTRINSIC
     for (int i1 = 0; i1 < n1; i1 += 8)
#else
     for (int i1 = 0; i1 < n1; ++i1)
#endif
     {
         for (int i2 = 0; i2 < n2; ++i2)
         // for (int i2 = 0; i2 < n2; i2+=2)
         {
#ifdef INTRINSIC
             size_t base_addrA = OFFSET_OF_INPUT(i1 - lowfil, i2 - lowfil, i3 - lowfil, n2, n1);
             // size_t base_addrB = OFFSET_OF_INPUT(i1 - lowfil, (i2+1) - lowfil, i3 - lowfil, n2, n3);

             SIMDFPTYPE x10A = _MM_LOAD(x1 + base_addrA);
             SIMDFPTYPE x20A = _MM_LOAD(x2 + base_addrA);
             SIMDFPTYPE tt1A = _MM_MUL(x10A, _MM_SET1(cy));
             SIMDFPTYPE tt2A = _MM_MUL(x20A, _MM_SET1(cy));

             SIMDFPTYPE coeff1, coeff2, x11, x22;
             size_t in_stride = IN_STRIDE1(n1);

             base_addrA += lowfil*in_stride;
             // base_addrB += lowfil*in_stride;

             size_t out_stride = OUT_STRIDE1(n1);
             size_t o_base_addrA = OFFSET_OF_OUTPUT(i1, i2, i3, n2, n1);
             // size_t o_base_addrB = OFFSET_OF_OUTPUT(i1, (i2+1), i3, n2, n3);

             _MM_PREFETCH1(x1 + base_addrA + (28 + 1)*in_stride);
             _MM_PREFETCH1(x2 + base_addrA + (28 + 1)*in_stride);
             _MM_PREFETCH1(y1 + o_base_addrA + 1*out_stride);
             _MM_PREFETCH1(y2 + o_base_addrA + 1*out_stride);

             BOOST_PP_REPEAT(29, MAD_FILTERA, 1);

             _mm512_storenrngo_pd(
               y1 + o_base_addrA, _MM_ADD(_MM_LOAD(y1 + o_base_addrA), tt1A));
             _mm512_storenrngo_pd(
               y2 + o_base_addrA, _MM_ADD(_MM_LOAD(y2 + o_base_addrA), tt2A));

             SIMDFPTYPE kstrt2_v = _MM_FMA(_MM_MUL(tt1A, x10A), tt2A, x20A);
             kstrt2 += _mm512_reduce_add_pd(kstrt2_v);

#else
             tt1=x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*cy;
             tt2=x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*cy;
             for (int l = lowfil; l <= lupfil; ++l)
             {
                 int j = (i2 + l);

                 tt1+=x1[OFFSET_OF_INPUT(i1 + (-lowfil),j + (-lowfil),i3 + (-lowfil), n2, n1)]*fil1[l + (-lowfil)][1] -
                      x2[OFFSET_OF_INPUT(i1 + (-lowfil),j + (-lowfil),i3 + (-lowfil), n2, n1)]*fil2[l + (-lowfil)][1];
                 tt2+=x2[OFFSET_OF_INPUT(i1 + (-lowfil),j + (-lowfil),i3 + (-lowfil), n2, n1)]*fil1[l + (-lowfil)][1]+
                      x1[OFFSET_OF_INPUT(i1 + (-lowfil),j + (-lowfil),i3 + (-lowfil), n2, n1)]*fil2[l + (-lowfil)][1];
             }
             y1[OFFSET_OF_OUTPUT(i1,i2,i3, n2, n1)]+=tt1;
             y2[OFFSET_OF_OUTPUT(i1,i2,i3, n2, n1)]+=tt2;

             kstrt2+=tt1*x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]+
                     tt2*x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)];
#endif
         }
      }
  }

  barrier(threadid);
  time_3 = __rdtsc();

  if (0 == threadid)
  {
  time_t1 += (time_2 - time_1);
  time_t2 += (time_3 - time_2);
  }
  time_1 = __rdtsc();

  // ! + (1/2) d^2/dz^2
  for (int i3 = i3_begin; i3 < i3_end; ++i3)
  {
      for (int i2 = 0; i2 < n2; ++i2)
      {
#ifdef INTRINSIC
          for (int i1 = 0; i1 < n1; i1 += 8)
#else
          for (int i1 = 0; i1 < n1; ++i1)
#endif
          {
#ifdef INTRINSIC
             size_t base_addr = OFFSET_OF_INPUT(i1 - lowfil, i2 - lowfil, i3 - lowfil, n2, n1);
             SIMDFPTYPE x10 = _MM_LOAD(x1 + base_addr);
             SIMDFPTYPE x20 = _MM_LOAD(x2 + base_addr);
             SIMDFPTYPE tt1 = _MM_MUL(x10, _MM_SET1(cx));
             SIMDFPTYPE tt2 = _MM_MUL(x20, _MM_SET1(cx));

             SIMDFPTYPE coeff1, coeff2, x11 = _mm512_undefined_pd(), x22 = _mm512_undefined_pd();
             size_t in_stride = VLEN;
             base_addr += lowfil;

             size_t out_stride = VLEN;
             size_t o_base_addr = OFFSET_OF_OUTPUT(i1, i2, i3, n2, n1);

             _MM_PREFETCH1(x1 + base_addr + (28 + L1_PREFETCH_DISTANCE)*in_stride);
             _MM_PREFETCH1(x2 + base_addr + (28 + L1_PREFETCH_DISTANCE)*in_stride);
             _MM_PREFETCH1(y1 + o_base_addr + L1_PREFETCH_DISTANCE*out_stride);
             _MM_PREFETCH1(y2 + o_base_addr + L1_PREFETCH_DISTANCE*out_stride);

#undef MAD_FILTER
#define MAD_FILTER(z, l, dummy) \
             coeff1 = _MM_SET1(fil1[0][l]), coeff2 = _MM_SET1(fil2[0][l]); \
             _MM_LOADU(x11, x1 + base_addr + l); \
             _MM_LOADU(x22, x2 + base_addr + l); \
             tt1 = _MM_MSUBR231(_MM_FMA(tt1, x11, coeff1), x22, coeff2); \
             tt2 = _MM_FMA(_MM_FMA(tt2, x22, coeff1), x11, coeff2);
             BOOST_PP_REPEAT(29, MAD_FILTER, 0);

             _mm512_storenrngo_pd(
               y1 + o_base_addr, _MM_ADD(_MM_LOAD(y1 + o_base_addr), tt1));
             _mm512_storenrngo_pd(
               y2 + o_base_addr, _MM_ADD(_MM_LOAD(y2 + o_base_addr), tt2));

             SIMDFPTYPE kstrt3_v = _MM_FMA(_MM_MUL(tt1, x10), tt2, x20);
             kstrt3 += _mm512_reduce_add_pd(kstrt3_v);
#else
               tt1=x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*cx;
               tt2=x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]*cx;

               for (int l = lowfil; l <= lupfil; ++l)
               {
                   int j = (i3 + l);
                
                   tt1+=x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),j + (-lowfil), n2, n1)]*fil1[l + (-lowfil)][2]-
                        x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),j + (-lowfil), n2, n1)]*fil2[l + (-lowfil)][2];
                   tt2+=x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),j + (-lowfil), n2, n1)]*fil1[l + (-lowfil)][2]+
                        x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),j + (-lowfil), n2, n1)]*fil2[l + (-lowfil)][2];
               } 
               y1[OFFSET_OF_OUTPUT(i1,i2,i3, n2, n1)]+=tt1;
               y2[OFFSET_OF_OUTPUT(i1,i2,i3, n2, n1)]+=tt2;

               kstrt3+=tt1*x1[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)]+
                       tt2*x2[OFFSET_OF_INPUT(i1 + (-lowfil),i2 + (-lowfil),i3 + (-lowfil), n2, n1)];
#endif
          }
      }
  }

  barrier(threadid);

  time_2 = __rdtsc();
  
  if (0 == threadid)
  time_t3 += (time_2 - time_1);

#pragma omp critical
  kstrten[0] += kstrt1;
  kstrten[1] += kstrt2;
  kstrten[2] += kstrt3;
} // END SUBROUTINE convolut_kinetic_per_T_k

#define _malloc_2M(X) \
  mmap(NULL, (X + 4096), PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED|MAP_POPULATE, -1, 0)

//   mmap(NULL, (X + 4096), PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED|MAP_HUGETLB|MAP_POPULATE, -1, 0)
void convolut_kinetic_per_T_k_wrapper(const unsigned int n1,
                                      const unsigned int n2,
                                      const unsigned int n3,
                                      const double *hgrid,
                                      const double *t3,
                                      const double *t4,
                                      double *y1,
                                      double *y2,
                                      double *kstrten,
                                      const double k1,
                                      const double k2,
                                      const double k3
                                      )
{
    // When we integrate these should be defined by the app already.
    const int lowfil = -14; 
    const int lupfil = 14;

    // In BigDFT x comes as a four dinemsional array. THe fourth dimension
    // only has 2 elements always. So it is split up as 2 3-d arrays
    static double *x1 = NULL;
    static double *x2 = NULL;

    if (x1 == NULL) 
    {
        x1 = (double *)_malloc_2M((IN_STRIDE3(n3) * IN_STRIDE21(n2, n1) + 8) * sizeof(double));
        x1 += 2;
    }

    if (x2 == NULL)
    {
        x2 = (double *)_malloc_2M((IN_STRIDE3(n3) * IN_STRIDE21(n2, n1) + 8)* sizeof(double));
        x2 += 2;
    }


#if 1
   unsigned long long int skirt_sum1 = 0; 
   unsigned long long int skirt_sum2 = 0; 
   unsigned long long int skirt_omp1 = 0; 
   // unsigned long long int skirt_omp2 = 0; 
   unsigned long long int skirt_time0; 
   unsigned long long int skirt_time1; 
   unsigned long long int skirt_time2; 
   // unsigned long long int skirt_time3; 
   // unsigned long long int skirt_time4; 
   // for (int ITERATION = 0; ITERATION < 1; ITERATION++)
   {
   skirt_time0 = __rdtsc();

#pragma omp parallel for
   for (int i3 = 0; i3 < n3; ++i3)
   {
      unsigned ID = omp_get_thread_num();
      if ((ID == 0) && (i3 == 0)) skirt_time1 = __rdtsc();

      double *t3_ptr = &(t3[i3*(n1*n2)]);
      double *t4_ptr = &(t4[i3*(n1*n2)]);
      for (int i2 = 0; i2 < n2; ++i2)
      {
         double *x1_ptr = &(x1[OFFSET_OF_INPUT(-lowfil, i2-lowfil, i3-lowfil, n2, n1)]);
         double *x2_ptr = &(x2[OFFSET_OF_INPUT(-lowfil, i2-lowfil, i3-lowfil, n2, n1)]);
         for (int i1 = 0; i1 < n1; i1+=16)
         {
            __m512 a0, a1, b0, b1;

            _mm_prefetch(&(t3_ptr[16]), _MM_HINT_T0);
            _mm_prefetch(&(t3_ptr[24]), _MM_HINT_T0);
            _mm_prefetch(&(t4_ptr[16]), _MM_HINT_T0);
            _mm_prefetch(&(t4_ptr[24]), _MM_HINT_T0);

            _mm_prefetch(&(x1_ptr[16]), _MM_HINT_ET0);
            _mm_prefetch(&(x1_ptr[24]), _MM_HINT_ET0);
            _mm_prefetch(&(x2_ptr[16]), _MM_HINT_ET0);
            _mm_prefetch(&(x2_ptr[24]), _MM_HINT_ET0);

            a0 = _mm512_extload_ps(&(t3_ptr[0]), _MM_UPCONV_PS_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            a1 = _mm512_extload_ps(&(t3_ptr[8]), _MM_UPCONV_PS_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            b0 = _mm512_extload_ps(&(t4_ptr[0]), _MM_UPCONV_PS_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            b1 = _mm512_extload_ps(&(t4_ptr[8]), _MM_UPCONV_PS_NONE, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            _mm512_extstore_ps(&(x1_ptr[0]), a0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
            _mm512_extstore_ps(&(x1_ptr[8]), a1, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
            _mm512_extstore_ps(&(x2_ptr[0]), b0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
            _mm512_extstore_ps(&(x2_ptr[8]), b1, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);

            t3_ptr += 16;
            t4_ptr += 16;
            x1_ptr += 16;
            x2_ptr += 16;
         }
      }
   }

   skirt_time2 = __rdtsc();


   // skirt_time4 = __rdtsc();

   skirt_sum1 += (skirt_time2 - skirt_time1);
   // skirt_sum2 += (skirt_time4 - skirt_time3);
   skirt_omp1 += (skirt_time1 - skirt_time0);
   // skirt_omp2 += (skirt_time3 - skirt_time2);

   // Case 1
    skirt_time1 = __rdtsc();
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 2
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 3
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 4
   #pragma omp parallel for
   for (int i3 = n3 + (-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 5
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 6
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 7
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 8
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = -lowfil; i1 < n1+(-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
               input_i1 = i1 + lowfil;

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 9
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 10
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 11
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 12
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 13
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 14
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 15
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 16
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 17
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = 0; i1 < (-lowfil); ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 <= ((-lowfil)-1))
               input_i1 = n1+i1+lowfil;

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 18
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 19
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 20
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = -lowfil; i2 < n2+(-lowfil); ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
               input_i2 = i2 + lowfil;

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 21
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 22
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 23
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = 0; i2 < (-lowfil); ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // else if (i2 <= ((-lowfil)-1))
               input_i2 = n2+i2+lowfil;

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 24
   #pragma omp parallel for
   for (int i3 = -lowfil; i3 < n3+(-lowfil); ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
               input_i3 = i3 + lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 25
   #pragma omp parallel for
   for (int i3 = 0; i3 < (-lowfil); ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // else if (i3 <= ((-lowfil)-1))
               input_i3 = n3+i3+lowfil;

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }

   // Case 26
   #pragma omp parallel for
   for (int i3 = n3+(-lowfil); i3 < n3+(-lowfil)+lupfil; ++i3)
   {
      for (int i2 = n2+(-lowfil); i2 < n2+(-lowfil)+lupfil; ++i2)
         for (int i1 = n1+(-lowfil); i1 < n1+(-lowfil)+lupfil; ++i1)
         {
            int input_i1, input_i2, input_i3;

            // else if (i1 >= (n1+(-lowfil)))
               input_i1 = i1 - (n1+(-lowfil));

            // else if (i2 >= (n2+(-lowfil)))
               input_i2 = i2 - (n2+(-lowfil));

            // else if (i3 >= (n3+(-lowfil)))
               input_i3 = i3 - (n3+(-lowfil));

            x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
            x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
               t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];
         }
   }
   skirt_time2 = __rdtsc();
   skirt_sum2 += (skirt_time2 - skirt_time1);
   }
   // printf("CPU freq %f Skirt Time (ms) = %f, %f, %f, %f\n", get_cpu_freq(), ((skirt_sum1)/get_cpu_freq())*1000, ((skirt_sum2)/get_cpu_freq())*1000, ((skirt_omp1)/get_cpu_freq())*1000, ((skirt_omp2)/get_cpu_freq())*1000);
#else
   long long skirt_start_time, skirt_end_time;
   skirt_start_time = __rdtsc();
// #pragma omp parallel for
   for (int i3 = 0; i3 < n3+(-lowfil)+lupfil; ++i3)
        for (int i2 = 0; i2 < n2+(-lowfil)+lupfil; ++i2)
            for (int i1 = 0; i1 < n1+(-lowfil)+lupfil; ++i1)
            {
                int input_i1, input_i2, input_i3;

                if (i1 <= ((-lowfil)-1))
                    input_i1 = n1+i1+lowfil;

                if ((i1 >= (-lowfil)) && (i1 < (n1+(-lowfil))))
                    input_i1 = i1 + lowfil;

                if (i1 >= (n1+(-lowfil)))
                    input_i1 = i1 - (n1+(-lowfil));

                if (i2 <= ((-lowfil)-1))
                    input_i2 = n2+i2+lowfil;

                if ((i2 >= (-lowfil)) && (i2 < (n2+(-lowfil))))
                    input_i2 = i2 + lowfil;

                if (i2 >= (n2+(-lowfil)))
                    input_i2 = i2 - (n2+(-lowfil));

                if (i3 <= ((-lowfil)-1))
                    input_i3 = n3+i3+lowfil;

                if ((i3 >= (-lowfil)) && (i3 < (n3+(-lowfil))))
                    input_i3 = i3 + lowfil;

                if (i3 >= (n3+(-lowfil)))
                    input_i3 = i3 - (n3+(-lowfil));

                // if ((i3 == (-lowfil)) && (i2 == (-lowfil)) && (i1 < (-lowfil)+10))
                   //  printf("%d %d %d %24.20f\n", input_i1, input_i2, input_i2, t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1]);

                x1[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
                      t3[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];

                x2[OFFSET_OF_INPUT(i1, i2, i3, n2, n1)] =
                      t4[(input_i3 * n2 * n1) + (input_i2*n1) + input_i1];

            }
     skirt_end_time = __rdtsc();
     printf("CPU freq %f Skirt Time = %f ms\n", get_cpu_freq(), ((skirt_end_time - skirt_start_time)/get_cpu_freq()) * 1000);
#endif

    //const int numtimes = 13 * 1008; //
   
    // It really does not make much sense to change it beyond 1
    const int numtimes = 1;
    nthreads = omp_get_max_threads();
    double clock_freq = get_cpu_freq();
    // printf("clock_freq = %f\n", clock_freq);

    unsigned long long int start_time, end_time, time_copy, time_convolute;
    unsigned long long int time_1, time_2;

#pragma omp parallel
    {
        unsigned threadid = omp_get_thread_num();
        for (int i = 0; i < numtimes; ++i) 
        {
            if (0 == threadid && 0 == i)
            {
              start_time = __rdtsc();
              time_init = time_t1 = time_t2 = time_t3 = 0; 
              time_copy = time_convolute = 0;
            }
            unsigned long long time_0 = __rdtsc();

            int i3_per_thread = (n3 + nthreads - 1)/nthreads;
            int i3_begin = i3_per_thread*threadid;
            int i3_end = MIN(i3_begin + i3_per_thread, n3);

            time_1 = __rdtsc();
            convolut_kinetic_per_T_k(n1, n2, n3, hgrid, x1, x2, y1, y2, kstrten, k1, k2, k3, threadid);
            time_2 = __rdtsc();
            if (0 == threadid)
            {
              time_copy += (time_1 - time_0);
              time_convolute += (time_2 - time_1);
            }
        }
    }
    end_time = __rdtsc();

    // time_to_convolute - total time for convolution
    // time_t1 - time for loop 1
    // time_t2 - time for loop 2
    // time_t3 - time for loop 3
    // skirt_omp1 - OMP overhead for skirt
    // skirt_sum1 - adding "core" skirt
    // skirt_sum2 - adding "peripheral" skirt
    printf("%f, %f, %f, %f, %f, %f, %f\n",
           ((time_convolute / numtimes)/clock_freq) * 1000,
           ((time_t1 / numtimes)/clock_freq) * 1000,
           ((time_t2 / numtimes)/clock_freq) * 1000,
           ((time_t3 / numtimes)/clock_freq) * 1000,
           ((skirt_omp1 / 1)/clock_freq) * 1000,
           ((skirt_sum1 / 1)/clock_freq) * 1000,
           ((skirt_sum2 / 1)/clock_freq) * 1000);
}
