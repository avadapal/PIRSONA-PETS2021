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

#include <unistd.h>
#include <fcntl.h>
#include <iostream>
#include <tuple>
//#include "common.h"
#include<cmath>
#include <algorithm>
#include<vector>
#include "dpf.h"
#include "aes.h"
#include "prg.h"
#include "sanitycheck.h"

template<size_t nitems> struct query { size_t user; dpf_key<nitems> dpfkey; };

int main(int argc, char ** argv)
{
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  // const size_t nitems = 1ULL << 10;
  // const size_t nqueries = 1ULL << 20; 
 
  ssize_t bytes;
  int fd = open("./queries", O_RDONLY);
  //bytes = read(fd, &nqueries, sizeof(size_t)); 
  size_t * q = new size_t[2 * nqueries];
  bytes = read(fd, q, 2 * nqueries * sizeof(size_t));
  close(fd);

  query<nitems> * q0 = new query<nitems>[nqueries];
  query<nitems> * q1 = new query<nitems>[nqueries];
  dpf_key<nitems> dpfkey[2] = { 0 }; 
  __m128i ** s0 = (__m128i**)malloc(nqueries * sizeof(__m128i*));
  uint8_t ** t0 = (uint8_t**)malloc(nqueries * sizeof(uint8_t*));

  __m128i ** s1 = (__m128i**)malloc(nqueries * sizeof(__m128i*));
  uint8_t ** t1 = (uint8_t**)malloc(nqueries * sizeof(uint8_t*));
 

  __m128i seeds[2];
  arc4random_buf(seeds, 2*sizeof(__m128i));
   std::vector<__m256i> rnd0_[4];
 std::vector<__m256i> rnd1_[4];
 auto [rnd0, rnd1] = gen_rands(aeskey, seeds, nitems);
   rnd0_[0] = rnd0; rnd0_[1] = rot(rnd0,1); rnd0_[2] = rot(rnd0,2); rnd0_[3] = rot(rnd0,3);
   rnd1_[0] = rnd1; rnd1_[1] = rot(rnd1,1); rnd1_[2] = rot(rnd1,2); rnd1_[3] = rot(rnd1,3); 
  for (size_t i = 0; i < nqueries; ++i)
  {
    gen(aeskey, q[2*i+1], dpfkey);

    // t0[i] = (uint8_t*)malloc(nitems * sizeof(uint8_t));
    // t1[i] = (uint8_t*)malloc(nitems * sizeof(uint8_t));
    //     if (posix_memalign((void**)&s0[i], sizeof(__m128i), nitems * sizeof(__m128i)))
    // {
    //  throw std::runtime_error("posix_memalign failed");
    // }

    //     if (posix_memalign((void**)&s1[i], sizeof(__m128i), nitems * sizeof(__m128i)))
    // {
    //  throw std::runtime_error("posix_memalign failed");
    // }


    // evalfull3(aeskey, dpfkey[1], s1[i], t1[i]); 
    // evalfull3(aeskey, dpfkey[0], s0[i], t0[i]);
    
    // const bool * buf_bool0 = reinterpret_cast<const bool*>(t0[i]);
    // const bool * buf_bool1 = reinterpret_cast<const bool*>(t1[i]);
    // const uint8_t * q0_ = reinterpret_cast<const uint8_t*>(buf_bool0); 
    // const uint8_t * q1_ = reinterpret_cast<const uint8_t*>(buf_bool1); 

    // uint64_t dot0 = dot_prod(rnd0_, 1, rnd1_, 1, q0_);
  
    // uint64_t dot1 = dot_prod(rnd0_, 1, rnd1_, 1, q1_);
    


    // for(size_t j = 0; j < nitems; ++j)
    // {
    //   std::cout << j << " ->: " << (int)t0[i][j] << " <> " << (int)t1[i][j] << std::endl;
    // }
 
    // std::cout << "dots = " << (dot0 - dot1) << " = " << dot0 << " - " << dot1 << std::endl;
    // #ifdef VERBOSE 
    // std::cout << "dpf: " << q[2*i] << " -- " << q[2*i + 1] << " -- " << dpfkey[0].leaf << " -- " << dpfkey[1].leaf << std::endl;
    // #endif

    q0[i] = { q[2*i], dpfkey[0] };
    q1[i] = { q[2*i], dpfkey[1] };
  }

  int fd0 = open("./queries.0", O_WRONLY, O_CREAT);
  //bytes = write(fd0, &nqueries, sizeof(size_t));
  bytes = write(fd0, q0, nqueries * sizeof(query<nitems>));
  close(fd0);

  int fd1 = open("./queries.1", O_WRONLY, O_CREAT);
  //bytes = write(fd1, &nqueries, sizeof(size_t));
  bytes = write(fd1, q1, nqueries * sizeof(query<nitems>));
  close(fd1);

  

  delete [] q;
  delete [] q0;
  delete [] q1;

  return 0;
}