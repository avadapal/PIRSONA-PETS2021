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
  
  int bad = 0;
  for (size_t i = 0; i < nqueries; ++i)
  {
    gen(aeskey, q[2*i+1], dpfkey);

    q0[i] = { q[2*i], dpfkey[0] };
    q1[i] = { q[2*i], dpfkey[1] };

    srand(i);
    size_t rr = rand() % 100; 

    if(rr < 90)
    {
       ++bad;
       arc4random_buf(&dpfkey[0].leaf, sizeof(profile<3*precision>));
       arc4random_buf(&dpfkey[1].leaf, sizeof(profile<3*precision>));
       q0[i] =   { q[i], dpfkey[0] };
       q1[i] =   { q[i], dpfkey[1] };
    }
  }

  std::cerr << "number of bad queries = " << bad << std::endl;
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
