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


#include <bsd/stdlib.h>
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
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(0, 1);

  {
   


  const  int rand_seed = (argc < 2) ? 1 : atoi(argv[1]);

  srand(rand_seed);  //(1);

  profile * users;
  posix_memalign((void**)&users, sizeof(__m256i), nusers * sizeof(profile));

  prof<precision> * usersA;
  posix_memalign((void**)&usersA, sizeof(__m256i), nusers * sizeof(prof<precision>));

  double p[DIM];
  for (size_t i = 0; i < nusers; ++i)
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
    users[i] = p;
     std:: cout << std::endl << " ----- " << std::endl;
    usersA[i] = p;
  }
  ssize_t bytes;
  int fd = open("./uprofilesA", O_WRONLY, O_CREAT);
  //bytes = write(fd, &nusers, sizeof(size_t));
  bytes = write(fd, users, nusers * sizeof(profile));


  ssize_t bytesA;
  int fdA = open("./uprofiles", O_WRONLY, O_CREAT);
  //bytes = write(fd, &nusers, sizeof(size_t));
  bytes = write(fdA, usersA, nusers * sizeof(prof<precision>));
  close(fdA);



    
    close(fd);
    free(users);


  prof<precision> * shares;
  posix_memalign((void**)&shares, sizeof(__m256i), nusers * sizeof(prof<precision>));
  
  prof<precision> * shares2;
  posix_memalign((void**)&shares2, sizeof(__m256i), nusers * sizeof(prof<precision>));
  

  //arc4random_buf(shares, nusers * sizeof(prof<precision>));
  


  for(size_t jj = 0; jj < nusers; ++jj){
    for(size_t j = 0; j < DIM; ++j) shares[jj].f64[j] = rand(); 
  }

  int fd0 = open("./uprofiles.0", O_WRONLY, O_CREAT);
  //bytes = write(fd0, &nusers, sizeof(size_t));
  bytes = write(fd0, shares, nusers * sizeof(prof<precision>));
  close(fd0);

  for (size_t i = 0; i < nusers; ++i)
  {
    shares[i] = usersA[i] - shares[i];
    shares2[i] = usersA[i] - shares[i];

    prof<precision> rec;
    for(size_t j = 0; j < DIM; ++j)
    {
      rec.f64[j] = shares[i].f64[j] + shares2[i].f64[j]; 
      std::cout << "rec: " << (double) rec.f64[j] << " ";
    } 
      std::cout << std::endl;
  }

  int fd1 = open("./uprofiles.1", O_WRONLY, O_CREAT);
  //bytes = write(fd1, &nusers, sizeof(size_t));
  bytes = write(fd1, shares, nusers * sizeof(prof<precision>));
  close(fd1);
  }

 




  // {

  // const size_t nusers = 1024;

  // prof<precision> * users;
  // posix_memalign((void**)&users, sizeof(__m256i), nusers * sizeof(prof<precision>));

  // double p[DIM];
  // for (size_t i = 0; i < nusers; ++i)
  // {
  //   double norm = 0.0L;
  //   for (size_t j = 0; j < DIM; ++j)  
  //   {
  //     p[j] = dist(e2);
  //     norm += p[j]*p[j];
  //   }
  //   norm = std::sqrt(norm);
  //   for (size_t j = 0; j < DIM; ++j)  
  //   {
  //     p[j] /= norm;
  //   }
  //   users[i] = p;


  // }

  // ssize_t bytes;
  // int fd = open("./uprofilesA", O_WRONLY, O_CREAT);
  // //bytes = write(fd, &nusers, sizeof(size_t));
  // bytes = write(fd, users, nusers * sizeof(prof<precision>));
  // close(fd);

  // prof<precision> * shares;
  // posix_memalign((void**)&shares, sizeof(__m256i), nusers * sizeof(prof<precision>));
  // arc4random_buf(shares, nusers * sizeof(prof<precision>));

  // int fd0 = open("./uprofiles.0", O_WRONLY, O_CREAT);
  // //bytes = write(fd0, &nusers, sizeof(size_t));
  // bytes = write(fd0, shares, nusers * sizeof(prof<precision>));
  // close(fd0);

  // for (size_t i = 0; i < nusers; ++i)
  // {
  //   shares[i] = users[i] - shares[i];
  // }

  // int fd1 = open("./uprofiles.1", O_WRONLY, O_CREAT);
  // //bytes = write(fd1, &nusers, sizeof(size_t));
  // bytes = write(fd1, shares, nusers * sizeof(prof<precision>));
  // close(fd1);

  // free(users);
  // free(shares);




  // }


 

  return 0;
}
