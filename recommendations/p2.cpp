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
#include <thread>
#include <iostream>
#include <chrono>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fstream>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/asio/ssl.hpp>
#include<vector>
 
#include<fstream>
#include "aes.h"
#include "network.h"
 

 #include "dpf.h"

 using namespace dpf;

typedef bool leaf_t;
typedef block<__m128i> node_t;
typedef AES_KEY prgkey_t;



using boost::asio::ip::tcp;
using boost::asio::ip::tcp;
namespace ssl = boost::asio::ssl;
typedef ssl::stream<tcp::socket> ssl_socket;
boost::asio::io_context io_context;
boost::asio::ssl::context ssl_context(boost::asio::ssl::context::tls); 
 prgkey_t prgkey;

   std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_read;
   std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_recv;
 //#include "isqrt.h"
 
 
 
  uint64_t key0 = 597349;
  uint64_t key1 = 121379;

 #include "compare.h"
 enum comparison_step
 {
  shifts_out0 = 0,
  shifts_out1,
  blinds_out0,
  blinds_out1,
  dpfs_out0,
  dpfs_out1,
  zero_blinds0_out,
  zero_blinds1_out,
  eval,
  randomness_in,
  sq_zero_randomness_in,
  send_sq_blinds_to_P1_out,
  send_sq_blinds_to_P0_out,
  comparison_num_steps,
};
size_t comparison_progress[comparison_step::comparison_num_steps] = { 0 };
//#include "p_other.h"



 
uint64_t share0[nrecords];
uint64_t share1[nrecords];
double rand0[nrecords];
double rand1[nrecords];
  
 
 void keep_polling(boost::asio::io_context& io_context)
 {
 

      while (
              comparison_progress[comparison_step::dpfs_out0]         < nrecords   ||
              comparison_progress[comparison_step::dpfs_out1]         < nrecords   || 
              comparison_progress[comparison_step::zero_blinds0_out]  < nrecords/4 ||
              comparison_progress[comparison_step::zero_blinds1_out]  < nrecords/4  
             )
    {
      io_context.reset();
      io_context.poll();
    }
 }




 
 
 block<__m256i> flags_blinds0[nrecords/256];  // Blinds for blinding the flag vectors
 block<__m256i> uv_blinds0[nrecords/4];       // Blinds for blinding the dot-product
 block<__m256i> flags_blinds1[nrecords/256];  // Blinds for blinding the flag vectors
 block<__m256i> uv_blinds1[nrecords/4];       // Blinds for blinding the dot-product

 uint64_t sq_blinds0[2][nrecords];
 uint64_t sq_blinds1[2][nrecords];

 uint64_t profiles_blinds0[2][nrecords + 1][DIM];
 uint64_t profiles_blinds1[2][nrecords + 1][DIM];
 
 uint64_t zero_alpha0[nrecords];
 uint64_t zero_alpha1[nrecords];

 block<__m256i> __zero_alpha0[nrecords/4];
 block<__m256i> __zero_alpha1[nrecords/4];


uint64_t blinds0[nrecords];
uint64_t blinds1[nrecords];


uint64_t profile_alpha0[nrecords];
uint64_t profile_alpha1[nrecords];

uint64_t alpha[nrecords]; 
 __m128i seed_profile_blinds0, seed_profile_blinds1;

 __m128i flags_seed[2], zeros_seed[2];


uint16_t shift0[nrecords];
uint16_t shift1[nrecords];
uint16_t shift[nrecords];
size_t new_nitems = 1ULL << 16;

  uint64_t dot_prod(uint64_t a[DIM], uint64_t b[DIM])
 {

   uint64_t dot_prod_ = 0;
   for(size_t dim = 0; dim < DIM; ++dim)
   {
      dot_prod_ += (a[dim] * b[dim]);
   }

   return dot_prod_;
 }


  void generate_blinds()
  {
    for(size_t j = 0; j < nrecords; ++j)
    {
      uint64_t gamma0, gamma1, alpha0, alpha1;
      arc4random_buf(&gamma0, sizeof(uint64_t));
      arc4random_buf(&gamma1, sizeof(uint64_t));
      arc4random_buf(&alpha0, sizeof(uint64_t));
      arc4random_buf(&alpha1, sizeof(uint64_t));

      blinds0[j] = gamma0 * gamma1 + alpha0;
      blinds1[j] = gamma0 * gamma1 + alpha1;

 
    }
  }



  void generate_profile_blinds()
  {
    for(size_t j = 0; j < nrecords; ++j)
    {
      profile_alpha0[j] = dot_prod(profiles_blinds0[0][j], profiles_blinds1[1][j]);
      profile_alpha1[j] = dot_prod(profiles_blinds0[1][j], profiles_blinds1[0][j]);
    } 
  }

 void generate_shifts()
 {

    for(size_t j = 0; j < nrecords; ++j)
    {
        shift0[j] = arc4random_uniform(new_nitems);
        shift1[j] = arc4random_uniform(new_nitems);
        shift[j]  = (shift0[j] + shift1[j]) % new_nitems;
        //auto alpha =  (shift0 + shift1) % new_nitems;
        //comparison_progress[comparison_step::shifts_gen] = j + 1;
    }
 }



 void generate_dpfs()
 {
   // size_t dpf_nitems = 1ULL << (intpart+fracpart-prefix_len+1);

   for(size_t j = 0; j < nrecords; ++j)
   {

    // while(comparison_progress[comparison_step::shifts_gen] < j + 1)
    // {
    //   std::this_thread::yield();
    // }
 
   auto [dpf0, dpf1] = dpf_key<>::gen(prgkey, new_nitems, shift[j]);

    dpf_from_P2_to_P0[j].dpfnitems = dpf0.nitems;
    //std::cout << "nitems = " << dpf_from_P2_to_P0[0].dpfnitems << std::endl;
    dpf_from_P2_to_P0[j].dpfdepth = dpf0.depth();
    for(size_t d = 0; d < dpf0.depth(); ++ d)
    {
      dpf_from_P2_to_P0[j].cw[d] = dpf0.cw[d];
      dpf_from_P2_to_P1[j].cw[d] = dpf1.cw[d];
    }
    dpf_from_P2_to_P0[j].root = dpf0.root;
    dpf_from_P2_to_P0[j].finalizer = dpf0.finalizer;

    dpf_from_P2_to_P1[j].dpfnitems = dpf1.nitems;
    dpf_from_P2_to_P1[j].dpfdepth = dpf1.depth();
    dpf_from_P2_to_P1[j].root = dpf1.root;
    dpf_from_P2_to_P1[j].finalizer = dpf1.finalizer;

    // const size_t to = dpf_nitems-1;
    //  auto alloc_size = dpf_key<>::interval_bytes(0, to);
    // __m256i  *query0 = reinterpret_cast<__m256i*>(std::aligned_alloc(alignof(__m256i), alloc_size));
    // __m256i  *query1 = reinterpret_cast<__m256i*>(std::aligned_alloc(alignof(__m256i), alloc_size));

    
    // dpf0.evalinterval(0, to, (bool*)query0, NULL, true);
    // dpf1.evalinterval(0, to, (bool*)query1, NULL, true);

    //      std::cout << "root0 = " << dpf0.root.mX[0] << " " << dpf0.root.mX[1] << std::endl;
    //      std::cout << "root1 = " << dpf1.root.mX[0] << " " << dpf1.root.mX[1] << std::endl;
    //      std::cout <<  query0[0][0] << " " << query0[0][1] << " " << query0[0][2] << " " << query0[0][3] << std::endl;
    //      std::cout <<  query1[0][0] << " " << query1[0][1] << " " << query1[0][2] << " " << query1[0][3] << std::endl;
         // // compute the parities
 

   // comparison_progress[comparison_step::dpfs_gen] = j + 1;
    }
   // std::cout << "dpfs generated\n";

 }
 


 
 void write_shifts_P0(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
 {
 
       
    async_write(sout, boost::asio::buffer(&shift0[j], sizeof(shift0[j])),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
           if(!ec)
           {
             if(j + 1 < nrecords) 
             {                  
              write_shifts_P0(j + 1, io_context, sout);
             }
           }          
           else
           {               
             write_shifts_P0(j , io_context, sout);
           }

          comparison_progress[comparison_step::shifts_out0] = j + 1;
        
        }); 
 }

 void write_shifts_P1(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
 {
    // while(comparison_progress[comparison_step::shifts_gen] < j + 1)
    // {
    //   std::this_thread::yield();
    // }

    async_write(sout, boost::asio::buffer(&shift1[j], sizeof(shift1[j])),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
           if(!ec)
           {
             if(j + 1 < nrecords) 
             {                  
              write_shifts_P1(j + 1, io_context, sout);
             }
           }          
           else
           {               
             write_shifts_P1(j , io_context, sout);
           }

          comparison_progress[comparison_step::shifts_out1] = j + 1;
        
        }); 
 } 

 
 
 


 void write_dpfs_to_P0(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
 {

   // std::cout << " write_dpfs_to_P0 " << std::endl;
 
     // while(comparison_progress[comparison_step::dpfs_gen] < j + 1) 
     //  {
     //    std::this_thread::yield();
     //  }

       //std::cout << "writing: " << dpf_from_P2_to_P0[j].dpfnitems << std::endl;
       async_write(sout, boost::asio::buffer(&dpf_from_P2_to_P0[j], sizeof(dpf_from_P2_to_P0[j])),
            [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
            
            {       
               if(!ec)
               {
                 if(j + 1 < nrecords) 
                 {                  
                  write_dpfs_to_P0(j + 1, io_context, sout);
                 }
               }          
               else
               {               
                 write_dpfs_to_P0(j , io_context, sout);
               }

              comparison_progress[comparison_step::dpfs_out0] = j + 1;
            
            }); 
 }


 void write_dpfs_to_P1(size_t j,  boost::asio::io_context& io_context, tcp::socket& sout)
 {

    //std::cout << " write_dpfs_to_P1 " << std::endl;
    
    // while(comparison_progress[comparison_step::dpfs_gen] < j + 1) 
    // { 
    //   std::this_thread::yield();
    // }
    
     //std::cout << "writing: " << dpf_from_P2_to_P1[j].dpfnitems << std::endl;
     async_write(sout, boost::asio::buffer(&dpf_from_P2_to_P1[j], sizeof(dpf_from_P2_to_P1[j])),
            [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
            
            {       
                if(!ec)
                {
                   if(j + 1 < nrecords) 
                   {
                    write_dpfs_to_P1(j + 1, io_context, sout);
                   }
                }          
                else
                {               
                  write_dpfs_to_P1(j , io_context, sout);
                }

              comparison_progress[comparison_step::dpfs_out1] = j + 1;
           
            }); 
 }



 void generate_randomness_for_squaring_and_zeroing()
 {
 
      for(size_t j = 0; j < nrecords; ++j)
      {
        zero_alpha0[j] = sq_blinds0[0][j] * sq_blinds1[1][j];
        zero_alpha1[j] = sq_blinds0[1][j] * sq_blinds1[0][j];

 
      }

      for(size_t j = 0; j < nrecords/4; ++j)
      {
         for(size_t i = 0; i < 4; ++i)
         {
          __zero_alpha0[j].mX[i] = zero_alpha0[j*4 + i];
          __zero_alpha1[j].mX[i] = zero_alpha1[j*4 + i];
         }
          
      }
 }
 
 



  void send_zero_blinds_to_P0(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
  {
  
 
 
     async_write(sout, boost::asio::buffer(&__zero_alpha0[j], sizeof(__zero_alpha0[j])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      {
                
          if(!ec)
          {
             if(j + 1 < nrecords/4) 
             {
               //std::cout << "send_zero_blinds_to_P0 " << j << std::endl; 
               send_zero_blinds_to_P0(j + 1, io_context, sout);
             }
          }          
          else
          {               
            send_zero_blinds_to_P0(j , io_context, sout);
          }
        
        comparison_progress[comparison_step::zero_blinds0_out] = j + 1;
      }); 
   }

  void send_zero_blinds_to_P1(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
  {
  
    
     
     async_write(sout, boost::asio::buffer(&__zero_alpha1[j], sizeof(__zero_alpha1[j])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      {       
          if(!ec)
          {
             
             if(j + 1 < nrecords/4) 
             {
              //std::cout << "send_zero_blinds_to_P1 " << j << std::endl;  
              send_zero_blinds_to_P1(j + 1, io_context, sout);
             }
    
          }          
          else
          {               
            send_zero_blinds_to_P1(j , io_context, sout);
          }

        comparison_progress[comparison_step::zero_blinds1_out] = j + 1;
     
      }); 
   }




  void send_profile_blinds_to_P0(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
  {
  
     
 
     async_write(sout, boost::asio::buffer(&profile_alpha0[j], sizeof(profile_alpha0[j])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      {
                
          if(!ec)
          {
             if(j + 1 < nrecords) 
             {
               //std::cout << "send_zero_blinds_to_P0 " << j << std::endl; 
               send_profile_blinds_to_P0(j + 1, io_context, sout);
             }
          }          
          else
          {               
            send_profile_blinds_to_P0(j , io_context, sout);
          }
        
        comparison_progress[comparison_step::zero_blinds0_out] = j + 1;
      }); 
   }

  void send_profile_blinds_to_P1(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
  {
  
 
     
     async_write(sout, boost::asio::buffer(&profile_alpha1[j], sizeof(profile_alpha1[j])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      {       
          if(!ec)
          {
             
             if(j + 1 < nrecords) 
             {
              //std::cout << "send_zero_blinds_to_P1 " << j << std::endl;  
              send_profile_blinds_to_P1(j + 1, io_context, sout);
             }
    
          }          
          else
          {               
            send_profile_blinds_to_P1(j , io_context, sout);
          }

        comparison_progress[comparison_step::zero_blinds1_out] = j + 1;
     
      }); 
   }

 


   void generate_randomness()
   {
     for(size_t j = 0; j < nrecords; ++j)
     {
        double r0 = ((double) rand() / (RAND_MAX));
        rand0[j] = r0;

        double r1 = ((double) rand() / (RAND_MAX));
        rand1[j] = r1;
        
        comparison_progress[comparison_step::randomness_in] = j + 1;
     }
   }


  void send_randomness_to_P0(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
  
       while(comparison_progress[comparison_step::randomness_in] < j + 1) 
       {
         std::this_thread::yield(); 
       }

       async_write(sout, boost::asio::buffer(&rand0[j], sizeof(rand0[j])),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
               if(j + 1 < nrecords) 
               {
                send_randomness_to_P0(j + 1, io_context, sout);
               }
      
            }          
            else
            {               
              send_randomness_to_P0(j , io_context, sout);
            }
        
        }); 
 
   }

  

    void send_randomness_to_P1(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
  
       while(comparison_progress[comparison_step::randomness_in] < j + 1) 
        {
          std::this_thread::yield(); 
        }
        
        async_write(sout, boost::asio::buffer(&rand1[j], sizeof(rand1[j])),
          [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {                
              if(j + 1 < nrecords) 
              {
               send_randomness_to_P1(j + 1, io_context, sout);
              }
            }          
            else
            {               
              send_randomness_to_P1(j , io_context, sout);
            }
        
        }); 
    }
 


 

 

int main(int argc, char* argv[])
{ 
  
  keys_read.reserve(nrecords);
  keys_recv.reserve(nrecords);
   //dpf_read.reserve(nprefixes + 1);
  AES_KEY key;
  for(size_t j = 0; j < nrecords; ++j)
  {
    arc4random_buf(&share0[j], sizeof(uint64_t));
    arc4random_buf(&share1[j], sizeof(uint64_t));
  }
    


   arc4random_buf(&flags_seed[0], sizeof(__m128i));
   arc4random_buf(&flags_seed[1], sizeof(__m128i));

   arc4random_buf(&zeros_seed[0], sizeof(__m128i));
   arc4random_buf(&zeros_seed[1], sizeof(__m128i));

   PRG(prgkey,  flags_seed[0], (__m128i *) flags_blinds0,  (nrecords/256 * sizeof(block<__m256i> )) / sizeof(__m128i));
   PRG(prgkey,  zeros_seed[0], (__m128i *) uv_blinds0,  (nrecords/4 * sizeof(block<__m256i> )) / sizeof(__m128i)); 
   PRG(prgkey,  flags_seed[1], (__m128i *) flags_blinds1, (nrecords/256 * sizeof(block<__m256i> )) / sizeof(__m128i));
   PRG(prgkey,  zeros_seed[1], (__m128i *) uv_blinds1, (nrecords/4 * sizeof(block<__m256i>)) / sizeof(__m128i)); 

 
   arc4random_buf(&seed_profile_blinds0 , sizeof(__m128i));
   arc4random_buf(&seed_profile_blinds1, sizeof(__m128i));
  
   for(size_t j = 0; j < nrecords/4; ++j)
   {
    sq_blinds0[0][4*j + 0] = uv_blinds0[j].mX[0];
    sq_blinds0[0][4*j + 1] = uv_blinds0[j].mX[1];
    sq_blinds0[0][4*j + 2] = uv_blinds0[j].mX[2];
    sq_blinds0[0][4*j + 3] = uv_blinds0[j].mX[3];

    sq_blinds0[1][4*j + 0] = uv_blinds1[j].mX[0];
    sq_blinds0[1][4*j + 1] = uv_blinds1[j].mX[1];
    sq_blinds0[1][4*j + 2] = uv_blinds1[j].mX[2];
    sq_blinds0[1][4*j + 3] = uv_blinds1[j].mX[3];
   }

   for(size_t i = 0; i < nrecords/256; ++i)
   {
     for(size_t j = 0; j < 256; ++j)
     {
       sq_blinds1[0][256*i + j] = (int)flags_blinds0[i].bits[j];
       sq_blinds1[1][256*i + j] = (int)flags_blinds1[i].bits[j];
     }
   }
 

   
   generate_randomness_for_squaring_and_zeroing();
   generate_shifts();
   generate_dpfs();
   generate_blinds();
   generate_profile_blinds();
   std::cout << "\n ------------------ Generation Done ------------------ \n";
 

 
  tcp::acceptor acceptor1(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1(acceptor1.accept());
  // std::cerr << "Listenting on port 2: " << PORT_P1_P2 << std::endl;


  tcp::acceptor acceptor1_a(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_a));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1_a(acceptor1_a.accept());
  // std::cerr << "Listenting on port 2: " << PORT_P1_P2_a << std::endl;

  tcp::acceptor acceptor1_b(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_b));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1_b(acceptor1_b.accept());
  // std::cerr << "Listenting on port 3: " << PORT_P1_P2_b << std::endl;

  tcp::acceptor acceptor1_c(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_c));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1_c(acceptor1_c.accept());
  // std::cerr << "Listenting on port 3: " << PORT_P1_P2_c << std::endl;
 

  tcp::acceptor acceptor1_d(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_d));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1_d(acceptor1_d.accept());
  // std::cerr << "Listenting on port 3: " << PORT_P1_P2_d << std::endl;

  tcp::acceptor acceptor1_e(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_e));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1_e(acceptor1_e.accept());
  // std::cerr << "Listenting on port 3: " << PORT_P1_P2_e << std::endl;

  tcp::acceptor acceptor1_f(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_f));
  // std::cout << " --- > " << std::endl;
  tcp::socket s1_f(acceptor1_f.accept());
  // std::cerr << "Listenting on port 3: " << PORT_P1_P2_f << std::endl;


  tcp::acceptor acceptor0(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0(acceptor0.accept());
  // std::cerr << "Listenting on port 5: " << PORT_P0_P2 << std::endl;

  tcp::acceptor acceptor0_a(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_a));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0_a(acceptor0_a.accept());
  // std::cerr << "Listenting on port 6: " << PORT_P0_P2_a << std::endl;

  tcp::acceptor acceptor0_b(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_b));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0_b(acceptor0_b.accept());
  // std::cerr << "Listenting on port 7: " << PORT_P0_P2_b << std::endl;

   tcp::acceptor acceptor0_c(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_c));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0_c(acceptor0_c.accept());
  // std::cerr << "Listenting on port 7: " << PORT_P0_P2_c << std::endl;

  tcp::acceptor acceptor0_d(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_d));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0_d(acceptor0_d.accept());
  // std::cerr << "Listenting on port 7: " << PORT_P0_P2_d << std::endl;
  

  tcp::acceptor acceptor0_e(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_e));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0_e(acceptor0_e.accept());
  // std::cerr << "Listenting on port 7: " << PORT_P0_P2_e << std::endl;

  tcp::acceptor acceptor0_f(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_f));
  // std::cout << " --- > " << std::endl;
  tcp::socket s0_f(acceptor0_f.accept());
  // std::cerr << "Listenting on port 7: " << PORT_P0_P2_f << std::endl;

 

   
   write(s0, boost::asio::buffer(&flags_seed[0], sizeof(__m128i)));
   write(s1, boost::asio::buffer(&flags_seed[1], sizeof(__m128i)));

   write(s0, boost::asio::buffer(&zeros_seed[0], sizeof(__m128i)));
   write(s1, boost::asio::buffer(&zeros_seed[1], sizeof(__m128i)));
 


   write(s0_f, boost::asio::buffer(&seed_profile_blinds0, sizeof(__m128i)));
   write(s1_f, boost::asio::buffer(&seed_profile_blinds1, sizeof(__m128i)));

   std::thread poller(keep_polling, std::ref(io_context)); 
   

   std::thread  zero_blinds_sender_P0(send_zero_blinds_to_P0, 0, std::ref(io_context),  std::ref(s0_e));  
   std::thread  zero_blinds_sender_P1(send_zero_blinds_to_P1, 0, std::ref(io_context),  std::ref(s1_e));
   
   std::thread  profile_blinds_sender_P0(send_profile_blinds_to_P0, 0, std::ref(io_context),  std::ref(s0_d));  
   std::thread  profile_blinds_sender_P1(send_profile_blinds_to_P1, 0, std::ref(io_context),  std::ref(s1_d));
 
   std::thread shifts_to_P0_writer(write_shifts_P0,0, std::ref(io_context), std::ref(s0_c));
   std::thread shifts_to_P1_writer(write_shifts_P1,0, std::ref(io_context), std::ref(s1_c));
 
   std::thread dpf_to_P0_writer(write_dpfs_to_P0, 0, std::ref(io_context), std::ref(s0_b));
   std::thread dpf_to_P1_writer(write_dpfs_to_P1, 0, std::ref(io_context), std::ref(s1_b));
 
  

 
   
    poller.join();
 
 
   zero_blinds_sender_P0.join();
   zero_blinds_sender_P1.join();

 
   
   profile_blinds_sender_P0.join();
   profile_blinds_sender_P1.join();
 

   shifts_to_P0_writer.join();
   shifts_to_P1_writer.join();

 

   dpf_to_P0_writer.join();
   dpf_to_P1_writer.join();
 
   
   

   return 0;
}

 