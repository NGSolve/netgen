#ifndef NETGEN_CORE_SIMD_HPP
#define NETGEN_CORE_SIMD_HPP

/**************************************************************************/
/* File:   simd.hpp                                                       */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include "ngcore_api.hpp"

#include "simd_generic.hpp"

#ifndef __CUDA_ARCH__

#ifdef NETGEN_ARCH_AMD64
#ifndef __SSE__
#define __SSE__
#endif
#include "simd_sse.hpp"
#endif

#ifdef __AVX__
#include "simd_avx.hpp"
#endif

#ifdef __AVX512F__
#include "simd_avx512.hpp"
#endif

#ifdef __aarch64__
#include "simd_arm64.hpp"
#endif

#endif // __CUDA_ARCH__

namespace ngcore
{
#ifndef __CUDA_ARCH__
#ifdef NETGEN_ARCH_AMD64
  /*
  NETGEN_INLINE auto HSum (SIMD<double,2> v1, SIMD<double,2> v2, SIMD<double,2> v3, SIMD<double,2> v4)
  {
    SIMD<double,2> hsum1 = my_mm_hadd_pd (v1.Data(), v2.Data());
    SIMD<double,2> hsum2 = my_mm_hadd_pd (v3.Data(), v4.Data());
    return SIMD<double,4> (hsum1, hsum2);
  }
  */
  
  NETGEN_INLINE auto GetMaskFromBits( unsigned int i )
  {
    return SIMD<mask64>::GetMaskFromBits(i);
  }
#endif
#endif // __CUDA_ARCH__
  
  NETGEN_INLINE void SIMDTranspose (SIMD<double,4> a1, SIMD<double,4> a2, SIMD <double,4> a3, SIMD<double,4> a4,
                                    SIMD<double,4> & b1, SIMD<double,4> & b2, SIMD<double,4> & b3, SIMD<double,4> & b4)
  {
    if constexpr (sizeof(a1.Lo()) == 16)
      {
        auto [h1,h2] = Unpack(a1,a2);
        auto [h3,h4] = Unpack(a3,a4);
        b1 = SIMD<double,4> (h1.Lo(), h3.Lo());
        b2 = SIMD<double,4> (h2.Lo(), h4.Lo());
        b3 = SIMD<double,4> (h1.Hi(), h3.Hi());
        b4 = SIMD<double,4> (h2.Hi(), h4.Hi());
      }
    else
      {
        b1 = SIMD<double,4> (a1[0], a2[0], a3[0], a4[0]);
        b2 = SIMD<double,4> (a1[1], a2[1], a3[1], a4[1]);
        b3 = SIMD<double,4> (a1[2], a2[2], a3[2], a4[2]);
        b4 = SIMD<double,4> (a1[3], a2[3], a3[3], a4[3]);
      }
  }
  
  template<int N>
  NETGEN_INLINE auto HSum (SIMD<double,N> s1, SIMD<double,N> s2)
  {
    return SIMD<double,2>(HSum(s1), HSum(s2));
  }

  template<int N>
  NETGEN_INLINE auto HSum (SIMD<double,N> s1, SIMD<double,N> s2, SIMD<double,N> s3, SIMD<double,N> s4 )
  {
    // return SIMD<double,4>(HSum(s1), HSum(s2), HSum(s3), HSum(s4));
    return SIMD<double,4>(HSum(s1, s2), HSum(s3,s4));
  }
}

#endif // NETGEN_CORE_SIMD_HPP
