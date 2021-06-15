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
#include <boost/bind.hpp>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>
 
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include "aes.h"
#include "network.h"
 

 #include "dpf.h"
 
 using namespace dpf;

typedef bool leaf_t;
typedef block<__m128i> node_t;
typedef AES_KEY prgkey_t;



using boost::asio::ip::tcp; 
using namespace std::chrono;

 


 #include "compare.h"
 
 enum comparison_step
 {
  profile_dots_done = 0,
  zero_computed,
  zero_blinded_in,
  zero_reconstructed,
  blinded_profiles_done,
  blinded_user_profiles_in,
  blinded_item_profiles_in,
  partial_sum_done,
  zero_blinded,
  shifts_in,
  shifts_to_Pb_out,
  shifted_diff_done,
  shifts_reconstruction_done,
  profile_alphas_in, 
  cary_blinds_in, 
  zero_blinded_0_in,
  zero_blinded_1_in,
  zero_blinded_0_out,
  zero_blinded_1_out,
  sq_blinded_0_in,
  sq_blinded_1_in,
  dpfs_in,
  extracted_in, 
  parities_computed,  
  comparison_num_steps, 
};

size_t comparison_progress[comparison_step::comparison_num_steps] = { 0 };

 
  double rand_from_user;

  __m128i seed_profiles; 
  __m128i flags_seed, zeros_seed;



  uint64_t scaled_rand;

  uint64_t carry_blinds_from_P2[nrecords];
 


 
 
  uint64_t partial_sum[nrecords];
  uint64_t blinded_partial_sum[nrecords];
  uint64_t blinded_b_recv[nrecords];
 
 
  uint64_t diff[nrecords];
  
  uint64_t uv[nrecords];
  uint64_t eval_vector[nrecords];
  block<__m256i> __eval_vector[nrecords/256]; 
  uint16_t shifts_from_P2[nrecords]; 
  uint16_t shifted_diff[nrecords];
  uint16_t shifts_reconstructed[nrecords];

  size_t new_nitems = 1ULL << 16;

  size_t shifted_input[nrecords];

  std::array<size_t, 2> shifted_bounds;
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  #include "dotproduct.h" 
  #include "zero-out.h"


/**
  Keeps resetting the io_context and calls the poll function, inorder to invoke the boost-asio handlers   

*/ 
void keep_polling(boost::asio::io_context& io_context)
{
    while(
           comparison_progress[comparison_step::zero_computed]               < nrecords ||       
           comparison_progress[comparison_step::extracted_in]                < nrecords ||
           comparison_progress[comparison_step::parities_computed]           < nrecords ||
           comparison_progress[comparison_step::shifts_reconstruction_done]  < nrecords ||
           comparison_progress[comparison_step::dpfs_in]                     < nrecords  
         )
  {
    io_context.reset();
    io_context.poll();
  }
}


/** 
 * Computes partial sums on the "profile dot product vector" vector
*/  

void compute_partial_sums()
{

   for(size_t i = 0; i < nrecords; ++i)
   {
     while(comparison_progress[comparison_step::zero_computed] < i + 1)
     {
      std::this_thread::yield();
     }

      uint64_t tmp_ = 0;
      
      for(size_t j = 0; j < i; ++j)
      {
        tmp_ += zeroed_vectors[j];
      }
      
      partial_sum[i] = tmp_;
      comparison_progress[comparison_step::partial_sum_done] = i + 1;
       
   } 
}


/**
  * Scales the uniformly selection random value \in \{0,1}, received from the user by last element of the partial sum array.
  * Calls the function "extract_bits" to extract the required bits from diff[j]

*/
void __extract_bits()
{

  for(size_t j = 0; j < nrecords; ++j)
  {
    while(comparison_progress[comparison_step::partial_sum_done] < nrecords)
     {
      std::this_thread::yield();
     }
     
     
     scaled_rand = rand_from_user * partial_sum[nrecords-1];
     diff[j]  =  extract_bits((partial_sum[j] - scaled_rand), 2*precision-fracpart, 2*precision+intpart-prefix_len);
     comparison_progress[comparison_step::extracted_in] = j + 1;
  }
}



 /** 
  * P0 and P1 read a DPFs at a random location \alpha
  * @param sin is the socket from which P0 and P1 receive the DPF 
 */
  void read_dpfs_from_P2(tcp::socket &sin)
  {
    for(size_t j = 0; j < nrecords; ++j)
    {
      read(sin, boost::asio::buffer(&dpf_read_from_P2[j], sizeof(dpf_read_from_P2[j])));
      comparison_progress[comparison_step::dpfs_in] = j + 1;
    }
  }

  
  /**
   * P0 and P1 receive the shares of the "random location", \alpha
   * @param sin is the socket from P0 and P1 receive the shares of the random location
  */
  void read_shifts_from_P2(tcp::socket &sin)
  {

    for(size_t j = 0; j < nrecords; ++j)
    {

       read(sin, boost::asio::buffer(&shifts_from_P2[j], sizeof(shifts_from_P2[j])));
       comparison_progress[comparison_step::shifts_in] = j + 1;
    }
  } 



  
  /*
   * Computes the shifted difference
   * If P0 and P1 are comparing a[j] and b[j], (diff[j] = a[j] - b[j])
   * "shifted diff" is the difference of diff[j] with the "random share" (i.e \alpha), received from P2
  */
  void compute_shifted_diffs()
  {
    
    for(size_t j = 0; j < nrecords; ++j)
    {
     
       while(
              comparison_progress[comparison_step::shifts_in]   < j + 1 || 
              comparison_progress[comparison_step::extracted_in] < j + 1
            )
       {
        std::this_thread::yield();
       }

      shifted_diff[j] = diff[j] - shifts_from_P2[j];
       
      comparison_progress[comparison_step::shifted_diff_done] = j + 1;
    }
  }



  /*
    sending shifts to do the reconstruct shifts.

    @param j
    @param io_context
    @param sout is the socket used to send the shifts
  */
  void send_shifts_to_Pb(size_t j, boost::asio::io_context& io_context, tcp::socket &sout)
  {

        while(
               comparison_progress[comparison_step::shifts_in] < j + 1
             )
       {
        std::this_thread::yield();
       }
      
      async_write(sout, boost::asio::buffer(&shifted_diff[j], sizeof(shifted_diff[j])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      {       
          if(!ec)
          {

             if(j + 1 < nrecords) 
             {             
              send_shifts_to_Pb(j + 1, io_context, sout);
             }
          }          
          else
          {               
            send_shifts_to_Pb(j , io_context, sout);
          }

           comparison_progress[comparison_step::shifts_to_Pb_out] = j + 1;
      
      }); 
  }


  /*
    Reconstructs the Shifts
    This computes (a - b - \alpha)
  */
  void reconstrutct_shifts(tcp::socket &sin)
  {
    for(size_t j = 0; j < nrecords; ++j)
    {
        read(sin, boost::asio::buffer(&shifts_reconstructed[j], sizeof(shifts_reconstructed[j])));
        shifts_reconstructed[j] = shifts_reconstructed[j] + shifted_diff[j];
        comparison_progress[comparison_step::shifts_reconstruction_done] = j + 1;
    }
  }

 


/**
 * Computes the required parities
 * 
 * "shifted_bounds" is the lookup table
 * 
 * Calls the "range_parities" function
 * Calls the "prefix_parities" function
 *
 * @param prgkey is the AES_KEY 
 * @param party indicates whether P0 or P1 is calling the function
*/
void  compute_parities(AES_KEY& prgkey, bool party)
{
  for(size_t j = 0; j < nrecords; ++j)
  {
    
      while(
          comparison_progress[comparison_step::dpfs_in]                    < j + 1 ||
          comparison_progress[comparison_step::shifts_reconstruction_done] < j + 1
         )
      {
        std::this_thread::yield();
      }

      shifted_bounds[0] =   (0 - shifts_reconstructed[j]); 
      shifted_bounds[1] =  (new_nitems/2 - shifts_reconstructed[j]); 

      shifted_bounds[0] = shifted_bounds[0] % new_nitems;
      shifted_bounds[1] = shifted_bounds[1] % new_nitems;

      size_t nitems_; 
      node_t root_;
      std::vector<node_t> cw_(dpf_read_from_P2[j].dpfdepth);
      std::array<node_t, nodes_per_leaf> finalizer_;
      
      nitems_ = dpf_read_from_P2[j].dpfnitems;//carry_dpf_read_from_P2[j].nitems; 
      root_ = dpf_read_from_P2[j].root; 
      for(size_t d = 0; d < dpf_read_from_P2[j].dpfdepth; ++d) cw_[d] = dpf_read_from_P2[j].cw[d];
     
      finalizer_ = dpf_read_from_P2[j].finalizer;
      dpf_key<> dpf_recv(std::move(nitems_),  std::move(root_), std::move(cw_), std::move(finalizer_), std::move(prgkey));

      auto bucket_parities  = range_parities(dpf_recv,  std::cbegin(shifted_bounds), std::cend(shifted_bounds));
      auto comparion_out    = get_poly_shares(bucket_parities, party, shifts_reconstructed[j]);

      comparison_progress[comparison_step::parities_computed] = j + 1;
  }
}

 
   
 

  


int main(int argc, char * argv[])
{  


  rand_from_user =  static_cast <float> (rand()) / static_cast <float> (RAND_MAX);


  uint64_t uprofile_reconstructed[DIM];
  uint64_t uprofile0[DIM], uprofile1[DIM];
  uint64_t iprofile0[DIM], iprofile1[DIM];
  
  uint64_t iprofile_reconstructed[nrecords][DIM];
  
  for(size_t dim = 0; dim < DIM; ++dim)
  {
    uprofile_reconstructed[dim] =  (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * scale;
    arc4random_buf(&uprofile0[dim], sizeof(uint64_t));
    uprofile1[dim] = uprofile_reconstructed[dim] - uprofile0[dim];

    for(size_t j = 0; j < nrecords; ++j)
    {
     iprofile_reconstructed[j][dim] = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * scale; 
    }
  }

 
  for(size_t dim = 0; dim < DIM; ++dim)
  {
    #if (PARTY == 0)
      uprofile[dim] = uprofile0[dim];
    #endif

    #if (PARTY == 1)
      uprofile[dim] = uprofile1[dim];
    #endif  

    
    for(size_t j = 0; j < nrecords; ++j)
    {   
      arc4random_buf(&iprofile0[dim], sizeof(uint64_t));
      iprofile1[dim] = iprofile_reconstructed[j][dim] - iprofile0[dim];
    
      #if (PARTY == 0)
        profiles[j][dim] = iprofile0[dim];
      #endif

      #if (PARTY == 1)
        profiles[j][dim] = iprofile1[dim];
      #endif  
    }
  }

  for(size_t dim = 0; dim < DIM; ++dim)
  {
    std::cout << double(uprofile_reconstructed[dim])/scale << " ";
  }

 std::cout << "\n\n ------------- " << std::endl;


 
  shifted_bounds[0] =  0; 
  shifted_bounds[1] =  new_nitems/2; 


   /*
    Generates some random eval vector for benchmarking purposes
   */
    for(size_t j = 0; j < nrecords; ++j)
    {
 
      uint8_t x;  
      arc4random_buf(&x, sizeof(uint8_t));
      x = x % 2;
      eval_vector[j] = 1;

      if(x == 0) eval_vector[j] = 0;

        #if(PARTY == 0)
         eval_vector[j] = 0;
        #endif
       __eval_vector[j/256].bits[j] = eval_vector[j];
     
     }

 
 

 bool party = true;

 #if (PARTY == 0)
    party = false;
 #endif

 

  const std::string host1 = (argc < 2) ? "127.0.0.1" : argv[1];
  const std::string host2 = (argc < 3) ? "127.0.0.1" : argv[2];
  const std::string host3 = (argc < 4) ? "127.0.0.1" : argv[3];
  

  AES_KEY prgkey;
 

  #include "connections.h"

  read(s2, boost::asio::buffer(&flags_seed, sizeof(__m128i)));
  read(s2, boost::asio::buffer(&zeros_seed, sizeof(__m128i)));



   
 // read(s2, boost::asio::buffer(&seed_zero, sizeof(__m128i)));
  read(s2_f, boost::asio::buffer(&seed_profiles, sizeof(__m128i)));
  

   PRG(prgkey, seed_profiles, (__m128i *) profile_blinds ,  (2 * DIM * (nrecords) * sizeof(uint64_t)) / sizeof(__m128i));
   PRG(prgkey, flags_seed, (__m128i *) flags_blinds ,  (nrecords/256 * sizeof(block<__m256i> )) / sizeof(__m128i));
   PRG(prgkey, zeros_seed, (__m128i *) uv_blinds ,  (nrecords/4 * sizeof(block<__m256i> )) / sizeof(__m128i)); 
 
   std::cout << "time measurement starts ......\n";

   std::thread poller(keep_polling, std::ref(io_context));
   
   std::thread profile_blinder(blind_profiles);
   std::thread user_profile_writer(write_blinded_user_profiles, 0, std::ref(io_context), std::ref(sb_f));
   std::thread user_profile_reader(read_blinded_user_profiles,std::ref(sb_f));
   std::thread item_profile_writer(write_blinded_item_profiles, 0, std::ref(io_context), std::ref(sb_g));
   std::thread item_profile_reader(read_blinded_item_profiles,std::ref(sb_g));
   std::thread profile_alpha_reader(read_profile_alphas_from_P2, std::ref(s2_d));
   std::thread dot_product_computer(compute_dot_products);

  

   std::thread zero_blinder(blinds_for_zeroing);  
   std::thread blinded_dot_prods_writer(write_blinded_dot_prods, 0, std::ref(io_context), std::ref(sb_a));
   std::thread blinded_flags_writer(write_blinded_flags, 0, std::ref(io_context), std::ref(sb_d)); 
   std::thread blinded_dot_prods_reader(read_blinded_dot_prods, std::ref(sb_a));
   std::thread blinded_flags_reader(read_blinded_flags, std::ref(sb_d));
   std::thread zero_blinds_from_P2_reader(read_zero_blinds_from_P2, std::ref(s2_e));
   std::thread zeroer(zeroing, std::ref(sb_b));

   std::thread partial_sum_computer(compute_partial_sums); 
   
   //DEBUGGING
   //std::thread zeroer_reconstructer(zero_reconstructing, std::ref(sb_b));
   
   std::thread dpf_reader(read_dpfs_from_P2, std::ref(s2_b));
   std::thread shifts_reader(read_shifts_from_P2, std::ref(s2_c));
   std::thread shifted_diff_computer(compute_shifted_diffs);
   std::thread shifts_sender_to_Pb(send_shifts_to_Pb, 0, std::ref(io_context), std::ref(sb_e));
   std::thread shifts_reconstructor(reconstrutct_shifts, std::ref(sb_e));
   std::thread extractor(__extract_bits);
   std::thread parity_computer(compute_parities, std::ref(prgkey), party);

 
     


    start = std::chrono::steady_clock::now();
    
    poller.join();

    profile_blinder.join();
    user_profile_writer.join();
    user_profile_reader.join();
    item_profile_writer.join();
    item_profile_reader.join();
    profile_alpha_reader.join();
    dot_product_computer.join();
    

    
    
    zero_blinder.join(); 
    blinded_dot_prods_writer.join();
    blinded_flags_writer.join();
    blinded_dot_prods_reader.join();
    blinded_flags_reader.join();
    zero_blinds_from_P2_reader.join();
    zeroer.join();
    

    partial_sum_computer.join();

    // DEBUGGING
    //zeroer_reconstructer.join(); 

    dpf_reader.join();
    shifts_reader.join();
    shifted_diff_computer.join();
    shifts_sender_to_Pb.join();
    shifts_reconstructor.join();
    extractor.join(); 
    parity_computer.join();

 
   std::cout << "evaluation is done ... \n"; 
   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
   std::cerr << "evaluation took "
              <<  double(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())/1000000
              << "seconds.\n";

  return 0;
}

 