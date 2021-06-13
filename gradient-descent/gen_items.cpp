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


#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <iostream>
#include <random>

#include "datatypesnp.h"
#include "fixed-np.h"
//#include "common.h"

int main(int argc, char ** argv)
{



{ 
  std::random_device rd;
  std::mt19937 e2(1);
  std::uniform_real_distribution<> dist(0, 1);
  
  // const size_t nqueries = 1ULL << 16;
  // const size_t nusers = 1024;
  // const size_t nitems = 1ULL << 10;
  const  int rand_seed = (argc < 2) ? 1 : atoi(argv[1]);
  srand(rand_seed);  //(1);

  profile * items;
  posix_memalign((void**)&items, sizeof(__m256i), nitems * sizeof(profile));

  prof<precision> * itemsA;
  posix_memalign((void**)&itemsA, sizeof(__m256i), nitems * sizeof(prof<precision>));

  double p[DIM];
  for (size_t i = 0; i < nitems; ++i)
  {
    double norm = 0.0L;
    for (size_t j = 0; j < DIM; ++j)	
    {
      //srand(j + i);
      p[j] = rand();  //dist(e2); //
      norm += p[j]*p[j];
    }
    norm = 1.0L / std::sqrt(norm);
    for (size_t j = 0; j < DIM; ++j)	
    {
      p[j] *= norm;
      std::cout << p[j] << " , ";
    }

     std:: cout << std::endl << " ----- " << std::endl;

    items[i]  = p;
    itemsA[i] = p;
  }


    for(size_t j = 0; j < nitems; ++j){
      for(size_t jj=0; jj < DIM; ++jj) std::cout << items[j].f64[jj] << " ";
      std::cout << std::endl;
    }



  int fd = open("./iprofilesA", O_WRONLY, O_CREAT);
  ssize_t bytes = write(fd, items, nitems * sizeof(profile));
  close(fd);

  free(items);


  int fdA = open("./iprofiles", O_WRONLY, O_CREAT);
  ssize_t bytesA = write(fdA, itemsA, nitems * sizeof(prof<precision>));
  close(fdA);

  free(itemsA);




}

  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(0, 1);

  {
  const  int rand_seed = (argc < 2) ? 1 : atoi(argv[1]);

  srand(rand_seed);  //(1);

  profile * items;
  posix_memalign((void**)&items, sizeof(__m256i), nitems * sizeof(profile));

  prof<precision> * itemsA;
  posix_memalign((void**)&itemsA, sizeof(__m256i), nitems * sizeof(prof<precision>));

  double p[DIM];
  for (size_t i = 0; i < nitems; ++i)
  {
    double norm = 0.0L;
    for (size_t j = 0; j < DIM; ++j)	
    {
       //srand(i+ j + 10000);
      p[j] =  rand(); 
      norm += p[j]*p[j];
    }
    norm = std::sqrt(norm);
    for (size_t j = 0; j < DIM; ++j)	
    {
      p[j] /= norm;
     std::cout << p[j] << " , ";
    }

    items[i] = p;
    std:: cout << std::endl << " ----- " << std::endl;
    itemsA[i] = p;
  }

  ssize_t bytes;
  int fd = open("./iprofilesA", O_WRONLY, O_CREAT);
  //bytes = write(fd, &nitems, sizeof(size_t));
  bytes = write(fd, items, nitems * sizeof(profile));


  ssize_t bytesA;
  int fdA = open("./iprofiles", O_WRONLY, O_CREAT);
  //bytes = write(fd, &nitems, sizeof(size_t));
  bytes = write(fdA, itemsA, nitems * sizeof(prof<precision>));
  close(fdA);



    
    close(fd);
    free(items);


  prof<precision> * shares;
  posix_memalign((void**)&shares, sizeof(__m256i), nitems * sizeof(prof<precision>));
  
  prof<precision> * shares2;
  posix_memalign((void**)&shares2, sizeof(__m256i), nitems * sizeof(prof<precision>));
  

  //arc4random_buf(shares, nitems * sizeof(prof<precision>));
  


  for(size_t jj = 0; jj < nitems; ++jj){
    for(size_t j = 0; j < DIM; ++j) shares[jj].f64[j] = rand(); 
  }

  int fd0 = open("./iprofiles.0", O_WRONLY, O_CREAT);
  //bytes = write(fd0, &nitems, sizeof(size_t));
  bytes = write(fd0, shares, nitems * sizeof(prof<precision>));
  close(fd0);

  for (size_t i = 0; i < nitems; ++i)
  {
    shares[i] = itemsA[i] - shares[i];
    shares2[i] = itemsA[i] - shares[i];

    prof<precision> rec;
    //std::cout << "SHARE:" << i << std::endl;
    for(size_t j = 0; j < DIM; ++j)
    {
      rec.f64[j] = shares[i].f64[j] + shares2[i].f64[j]; 
    //  std::cout << (double)shares[i].f64[j] << " + " <<  (double)shares2[i].f64[j] << " " << std::endl ;
     // std::cout << "rec: " << (double) rec.f64[j] << " ";
    } 
      std::cout << std::endl << std::endl;
  }

  int fd1 = open("./iprofiles.1", O_WRONLY, O_CREAT);
  //bytes = write(fd1, &nitems, sizeof(size_t));
  bytes = write(fd1, shares, nitems * sizeof(prof<precision>));
  close(fd1);
  }

  return 0;
}
