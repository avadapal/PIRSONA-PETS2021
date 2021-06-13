/*
PIRSONA: A 1-private 4-party gradient descent algorithm to provide recommendations

Copyright (C) 2021  Adithya Vadapalli, Ryan Henry

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DPF_H__
#define DPF_H__

#include <omp.h>
#include <cstdint>
#include <cmath>
#include <x86intrin.h>
#include <bsd/stdlib.h>

#include "aes.h"
#include "datatypes.h"
#include "prg.h"

#define L 0
#define R 1



static const __m128i lsb_mask[3] = {
  _mm_setzero_si128(),                                     // 0b00...0000
  _mm_set_epi64x(0,1)                                      // 0b00...0001
};
static const __m128i lsb_mask_inv = _mm_set_epi64x(-1,-2); // 0b11...1110
static const __m128i if_mask[2] = {
  _mm_setzero_si128(),                                     // 0b00...0000
  _mm_set1_epi8(-1)                                        // 0b11...1111
};

inline uint8_t _mm_getlsb_si128(const __m128i & x)
{
  __m128i vcmp = _mm_xor_si128(_mm_and_si128(x, lsb_mask[1]), lsb_mask[1]);
  return static_cast<uint8_t>(_mm_testz_si128(vcmp, vcmp));
}
inline __m128i _mm_clearlsb_si128(const __m128i & x)
{
  return _mm_and_si128(x, lsb_mask_inv);
}
inline __m128i _mm_setlsb_si128(const __m128i & x, bool b = true)
{
  return _mm_or_si128(_mm_clearlsb_si128(x), lsb_mask[b ? 1 : 0]);
}
inline __m128i _mm_xorif_si128(const __m128i & x, const __m128i & y, bool b)
{
  return _mm_xor_si128(x, _mm_and_si128(y, if_mask[b ? 1 : 0]));
}

struct cwbits { private: uint8_t t[2]; public: uint8_t & operator[](bool b) { return t[b ? 1 : 0]; } };

template<size_t nitems> struct dpf_key
{
  static_assert(nitems % 64 == 0, "nitems must be a multiple of 64");
  static constexpr size_t depth = static_cast<size_t>(std::ceil(std::log2(nitems)));

  __m128i root;
  __m128i cw[depth];
  cwbits t[depth];
  profile<3*precision> leaf;
};

inline void expand(AES_KEY & aeskey, const __m128i & seed, __m128i s[2],
  uint8_t t[2])
{
  const __m128i seedL = _mm_clearlsb_si128(seed);
  const __m128i seedR = _mm_setlsb_si128(seed);

  s[L] = seedL;
  s[R] = seedR;

  AES_ecb_encrypt_blks(s, 2, &aeskey); 

  s[L] = _mm_xor_si128(s[L], seedL);
  t[L] = _mm_getlsb_si128(s[L]);
  s[L] = _mm_clearlsb_si128(s[L]);

  s[R] = _mm_xor_si128(s[R], seedR);
  t[R] = _mm_getlsb_si128(s[R]);
  s[R] = _mm_clearlsb_si128(s[R]);
}

inline void nextlayer(AES_KEY & aeskey, __m128i s[2], uint8_t t[2], uint8_t bit,
  __m128i & cw, cwbits & cwt)
{
  __m128i s0[2], s1[2];
  uint8_t t0[2], t1[2];

  expand(aeskey, s[0], s0, t0);
  expand(aeskey, s[1], s1, t1);

  const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
  cw = _mm_clearlsb_si128(_mm_xor_si128(s0[lose], s1[lose]));
  cwt[L] = t0[L] ^ t1[L] ^ bit ^ 1;
  cwt[R] = t0[R] ^ t1[R] ^ bit;

  s[L] = _mm_xorif_si128(s0[keep], cw, t[L]);
  t[L] = t0[keep] ^ (t[L] & cwt[keep]);
  s[R] = _mm_xorif_si128(s1[keep], cw, t[R]);
  t[R] = t1[keep] ^ (t[R] & cwt[keep]);
}

inline void leaflayer(AES_KEY & aeskey, __m128i s[2], uint8_t t[2],
  profile<3*precision> & lastL, profile<3*precision> & lastR)
{
  constexpr size_t len = sizeof(profile<precision>) / sizeof(__m128i);

  profile<3*precision> s0, s1;
  PRG(aeskey, s[L], (__m128i*)&s0, dim128);
  PRG(aeskey, s[R], (__m128i*)&s1, dim128);
 

 
   // srand(1);
   // for(size_t j = 0; j < DIM; ++j) 
   //  {
   //    uint64_t x;
   //    x= 10;
   //    lastL.u64[j] = x;//rand(); 
   //  }
  arc4random_buf(&lastL, sizeof(profile<precision>));  
  
 lastR = (t[L] ? (s1 - s0) - lastL :  (s0 - s1) - lastL);

 lastL = -lastL; lastR = -lastR;

  #ifdef VERBOSE
 std::cout << "lastL = " << lastL << std::endl;
 std::cout << "lastR = " << lastR << std::endl; 
 std::cout << "lastR + lastL = " <<  lastR + lastL << std::endl;
 #endif

}

 

template <size_t nitems>
void gen(AES_KEY & aeskey, size_t point, dpf_key<nitems> dpfkey[2])
{
  __m128i s[2], s0[2], s1[2];
  uint8_t t[2], t0[2], t1[2];
 


     srand(1);
     s[0][0] = rand(); s[0][1] = rand();
     s[1][0] = rand(); s[1][1] = rand();
    //arc4random_buf(s, 2 * sizeof(__m128i));
 

  t[0] = _mm_getlsb_si128(s[0]);
  dpfkey[0].root = s[0];
  t[1] = !t[0];
  dpfkey[1].root = _mm_setlsb_si128(s[1], t[1]);

  for (size_t i = 0; i < dpf_key<nitems>::depth; ++i)
  {
    const uint8_t bit = (point >> (dpf_key<nitems>::depth - i - 1)) & 1U;
    nextlayer(aeskey, s, t, bit, dpfkey[0].cw[i], dpfkey[0].t[i]);
  }
  
  leaflayer(aeskey, s, t, dpfkey[0].leaf, dpfkey[1].leaf);
 

  memcpy(&dpfkey[1].cw, &dpfkey[0].cw, sizeof(dpf_key<nitems>::cw));
  memcpy(&dpfkey[1].t, &dpfkey[0].t, sizeof(dpf_key<nitems>::t));
}

template <size_t nitems>
inline void dolayer(AES_KEY & aeskey,__m128i * s, uint8_t * t,
  const size_t stepsize, dpf_key<nitems> dpfkey, size_t i, __m128i child[2],
  uint8_t ts[2])
{
  for (size_t j = 0; j < nitems; j += 2*stepsize)
  {
    expand(aeskey, s[j], child, ts);
    if (j+stepsize < nitems)
    {
      s[j+stepsize] = _mm_xorif_si128(child[R], dpfkey.cw[i], t[j]);
      t[j+stepsize] = ts[R] ^ dpfkey.t[i][R] & t[j];
    }
    s[j] = _mm_xorif_si128(child[L], dpfkey.cw[i], t[j]);
    t[j] = ts[L] ^ dpfkey.t[i][L] & t[j];
  }
}

 

/**
 * Does an Evalfull of DPF
 * @param aeskey is the PRGKEY
 * @param dpfkey is the DPFKEY
 * @param *s is the output of evalfull
*/ 
template <size_t nitems>
inline void evalfull3(AES_KEY & aeskey, dpf_key<nitems> & dpfkey,
__m128i * s, uint8_t * t)
{
  constexpr size_t depth = dpf_key<nitems>::depth;

  __m128i * s_[2] = { s, s + nitems/2 };
  uint8_t * t_[2] = { t, t + nitems/2 };

  int curlayer = depth % 2;

  s_[curlayer][0] = dpfkey.root;
  t_[curlayer][0] = _mm_getlsb_si128(dpfkey.root);

  __m128i child[2];
  uint8_t ts[2];
  for (size_t i = 0; i < depth; ++i)
  {
    curlayer = 1 - curlayer;
    const size_t itemnumber = std::max((nitems / (1ULL << (depth - i))), 1ULL);
    for (size_t j = 0; j < itemnumber; ++j)
    {
      expand(aeskey, s_[1-curlayer][j], child, ts);
      s_[curlayer][2*j] = _mm_xorif_si128(child[L], dpfkey.cw[i], t_[1-curlayer][j]);
      t_[curlayer][2*j] = ts[L] ^ dpfkey.t[i][L] & t_[1-curlayer][j];
      if (2*j+1 < 2*itemnumber)
      {
        s_[curlayer][2*j+1] = _mm_xorif_si128(child[R], dpfkey.cw[i], t_[1-curlayer][j]);
        t_[curlayer][2*j+1] = ts[R] ^ dpfkey.t[i][R] & t_[1-curlayer][j];
      }
    }
  }
}




template <size_t nitems, bool negate = false>
inline void dpf_mux(AES_KEY & aeskey,
  dpf_key<nitems> & dpfkey, profile<precision> * iprofiles,
  profile<precision> & itm, __m128i * s, uint8_t * t, fixed_t<2*precision>& cnt2)
{
  //evalfull(aeskey, dpfkey, s[0], t[0]); uint8_t * tt = t[0];
  // std::cout << "mux called" << std::endl;
  evalfull3(aeskey, dpfkey, s, t); 
  uint8_t * tt = t;// t[dpf_key<nitems>::depth % 2];
  memset(&itm, 0, sizeof(profile<precision>));
  uint64_t cnt = 0;
  
  for (size_t i = 0; i < nitems; ++i)
  {
    if (tt[i])
    {
      if constexpr(negate)
      {
        itm -= iprofiles[i];
      }
      else
      {
        itm += iprofiles[i];
      }
      cnt++;
    }
  }
  cnt2 = negate ? -static_cast<double>(cnt) : static_cast<double>(cnt);
    //return cnt;
}



/**
 * DPF based multiplexer
 * @param aeskey is the PRG key
 * @param dpfkey is the DPF key
 * @param *iprofiles is the array of profiles
 * @param itm is the output of the DPF based multiplexer
 * @param *t is the flag outputs the DPF evalfull
 * @param cnt2 is tha parities

*/
template <size_t nitems, bool negate = false>
inline void dpf_mux2(AES_KEY & aeskey,
  dpf_key<nitems> & dpfkey, profile<precision> * iprofiles,
  profile<precision> & itm, __m128i * s, uint8_t * t, fixed_t<2*precision>& cnt2, size_t j)
{
  
 
  uint8_t * tt = t;// t[dpf_key<nitems>::depth % 2];
  
  memset(&itm, 0, sizeof(profile<precision>));
  
  uint64_t cnt = 0;
  
  for (size_t i = 0; i < nitems; ++i)
  {  
   if (tt[i])
    {
      if constexpr(negate)
      {
        itm -= iprofiles[i];
      }
      else
      {
        itm += iprofiles[i];
      } 
      
      cnt++;
    }
  }
   cnt2 = negate ? -static_cast<double>(cnt) : static_cast<double>(cnt);
}



template <size_t nitems, bool negate = false>
inline void dpf_demux(AES_KEY & aeskey, dpf_key<nitems> & dpfkey,
  const profile<3*precision> & finalblock, profile<3*precision> * iprofiles,
  const __m128i * s, const uint8_t * t)
{

     profile<3*precision> tmp;
     //std::cout << "dpf_demux: " << std::endl;
     for (size_t i = 0; i < nitems; ++i)
     {
       PRG(aeskey, s[i], (__m128i *)&tmp,  dim128);
       
       if (t[i])
       {
        tmp -= finalblock;
       }
       else{

        //std::cout << i << ": " << tmp << std::endl;
       }

       if constexpr(negate) { iprofiles[i] = iprofiles[i] + tmp; /*std::cout << "N -> " << iprofiles[i] << std::endl;*/ }
       else { iprofiles[i] = iprofiles[i] - tmp; /*std::cout << "P-> " << iprofiles[i] << std::endl;*/  }      
     }
}



/**
 * DPF based multiplexer
 * @param dpfkey a DPF key
 * @param *tmp
 * @param finalblock 
 * @param *iprofiles is the array of profiles

*/
template <size_t nitems, bool negate = false>
inline void dpf_demux(dpf_key<nitems> & dpfkey,  profile<3*precision> * tmp,
  const profile<3*precision> & finalblock, profile<3*precision> * iprofiles,
  const __m128i * s, const uint8_t * t)
{

     for (size_t i = 0; i < nitems; ++i)
     { 
       if (t[i])
       {
        tmp[i] -= finalblock;
       }
       else{

        //std::cout << i << ": " << tmp << std::endl;
       }

       if constexpr(negate) { iprofiles[i] = iprofiles[i] + tmp[i]; /*std::cout << "N -> " << iprofiles[i] << std::endl;*/ }
       else { iprofiles[i] = iprofiles[i] - tmp[i]; /*std::cout << "P-> " << iprofiles[i] << std::endl;*/  }      
     }
}

#endif