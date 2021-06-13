#include <thread>
#include <iostream>
#include <chrono>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fstream>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>
#include<vector>
#include <cmath>
#include<fstream>
#include "aes.h"
#include "network.h"
 

 #include "dpf2.h"

 using namespace dpf;

typedef bool leaf_t;
typedef block<__m128i> node_t;
typedef AES_KEY prgkey_t;



using boost::asio::ip::tcp;
using namespace std::chrono;


 
 std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_recv;
 std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_b;
 std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_other;


 #include "isqrt.h"
 profile<precision> * uprofiles;
 enum norm_step
 {
  norm_bvrs_in = 0,
  dotprod_bvrs_in,
  profile_adjust,
  prenorm_profiles_in,
  prenorm_profiles_gen,
  unorms_out,
  msb_in,
  blinded_msb_in,
  blinded_msb_gen,
  dpf_gen,
  dpfs_in,
  s_other_in,
  s_compute,
  z_compute,
  z_gen,
  zz_in,
  z_added_done,
  blind_cnt_gen,
  blind_cnt_in,
  blind_sshare_gen,
  blind_sshare_in,
  sshare_gen,
  norm_done,
  msb_prod_gen,
  offset_in,
  offset_subtracted,
  normalized_out,
  normalized_in,
  norm_num_steps,
};

size_t norm_progress[norm_step::norm_num_steps] = { 0 };

 


  uint64_t key0 = 597349;
  uint64_t key1 = 121379;

 
int64_t cnt[nprofiles] = {0};
uint64_t s[nprofiles];
uint64_t s_other[nprofiles];
uint64_t blind_cnt[nprofiles];
uint64_t blind_sshare[nprofiles];
uint64_t blinded_cnt_recv[nprofiles]; 
uint64_t blinded_sshare_recv[nprofiles]; 
int64_t z[nprofiles];
int64_t zz[nprofiles];
uint64_t s_s_other[nprofiles];
uint64_t sshare0[nprofiles];
uint64_t sshare1[nprofiles];
uint32_t target[nprofiles];
int64_t z_added[nprofiles];
int64_t z_recv[nprofiles];


uint64_t profile_msb[nprofiles][DIM];
uint64_t profile_msb_recv[nprofiles][DIM];
uint64_t profile_msb_this[nprofiles][DIM];
uint64_t profile_msb_written[nprofiles][DIM];
uint64_t blinded_profile_msb[nprofiles][DIM];
uint64_t blinded_profile_msb_recv[nprofiles][DIM];
uint64_t profile_prod[DIM][nprofiles];
uint64_t profile_read[DIM][nprofiles];


uint64_t unorm_msb[nprofiles];
uint64_t unorm_msb_recv[nprofiles];
uint64_t unorm_msb_this[nprofiles];
uint64_t unorm_msb_written[nprofiles];
uint64_t blinded_unorm_msb[nprofiles];
uint64_t blinded_unorm_msb_recv[nprofiles];
uint64_t msb_prod[nprofiles];
uint64_t msb_read[nprofiles];

uint64_t reduced[nprofiles];
uint128_t _128_reduced[nprofiles];

uint64_t * query_recv[nprofiles]; 

profile<3 * precision> * normalized_profile_this  = (profile<3 * precision> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile<3 * precision>));
profile<3 * precision> * normalized_profile_other = (profile<3 * precision> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile<3 * precision>));
 
norm_bvrs * norm_bvrs_recv = (norm_bvrs *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(bvrs));

profile<param> * blinded_prenorm_profile = (profile<param> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile<param>));
profile<param> *  prenorm_uprofiles_;
profile<param> *  prenorm_uprofiles = (profile<param> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile<param>));
profile<param> *  blinded_prenorm_profile_out = (profile<param> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile<param>));

profile128<param> *  _128_blinded_prenorm_profile = (profile128<param> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile128<param>));
profile128<param> *  _128_prenorm_uprofiles = (profile128<param> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile128<param>));
profile128<param> *  _128_blinded_prenorm_profile_out = (profile128<param> *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(profile128<param>));

profile<3 * precision> * prenorm_uprofiles_3P;
fixed_t<2 * param> * unorm;
fixed_t<2 * param, uint128_t> * _128_unorm;
 void keep_polling(boost::asio::io_context& io_context)
{

   while (norm_progress[norm_step::msb_prod_gen] < nprofiles || norm_progress[norm_step::offset_subtracted] < nprofiles || 
        norm_progress[norm_step::dpf_gen] < nprofiles || norm_progress[norm_step::dpfs_in] < nprofiles || norm_progress[norm_step::s_compute] < nprofiles
        || norm_progress[norm_step::s_other_in] < nprofiles || norm_progress[norm_step::sshare_gen] < nprofiles ||  norm_progress[norm_step::blind_sshare_gen] < nprofiles
        || norm_progress[norm_step::blind_cnt_gen] < nprofiles || norm_progress[norm_step::blind_sshare_in] < nprofiles ||
           norm_progress[norm_step::blind_cnt_in] < nprofiles
        )
  {
    io_context.reset();
    io_context.poll();
  }
}


/**
 * @param io_context
 * @param sin is the socket used to read prenormalized profiles
*/
void read_prenorm_profiles(boost::asio::io_context & io_context, tcp::socket & sin)
{ 
  for (size_t i = 0; i < nprofiles; ++i)
  {
    read(sin, boost::asio::buffer(&_128_blinded_prenorm_profile[i], sizeof(profile128<param>)));
    norm_progress[norm_step::prenorm_profiles_in] = i + 1;    
  }
}



/**
 * blinds prenormalized profiles
*/
void blind_prenorm_profiles()
{ 

  if (posix_memalign((void**)&blinded_prenorm_profile_out, sizeof(__m256i),
    nprofiles * sizeof(profile<param>)))
  {
    throw std::runtime_error("Failed allocating blinded_profile");
  }

  for (size_t i = 0; i < nprofiles; ++i)
  {

    while (norm_progress[norm_step::dotprod_bvrs_in] < i + 1) 
     {
 
       std::this_thread::yield(); 
     }

    blinded_prenorm_profile_out[i]      = prenorm_uprofiles[i] + dotprod_recv[i].u_blind;
    _128_blinded_prenorm_profile_out[i] = _128_prenorm_uprofiles[i] + dotprod_recv[i]._128_u_blind;

    norm_progress[norm_step::prenorm_profiles_gen] = i + 1;
 
  }

}


/**
 * writes the prenormalized blinded profiles
 * @param i
 * @param io_context
 * @param sout is the socket used to write the pre-normalized blinded profiles
*/
void write_prenorm_profiles(size_t i, boost::asio::io_context & io_context, tcp::socket & sout)
{ 

  while(norm_progress[norm_step::prenorm_profiles_gen] < i + 1)
  { 
 
    std::this_thread::yield();
  }

   async_write(sout, boost::asio::buffer(&_128_blinded_prenorm_profile_out[i], sizeof(profile128<param>)),
        [i, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        { 
            if(!ec)
            {
               if(i + 1 < nprofiles){ 
                write_prenorm_profiles(i + 1, io_context, sout);
              }
            }          
            else
            {
              write_prenorm_profiles(i , io_context, sout);
            }
        
        }); 
}

 


/**
 * @param sin is the socket used to receive blinded dotproducts
*/
void read_dotprod_bvrs(tcp::socket& sin)
{
  for(size_t j = 0; j < nprofiles; ++j)
  {
 
    read(sin, boost::asio::buffer(&dotprod_recv[j], dotprod_stuff::size));

    norm_progress[norm_step::dotprod_bvrs_in] = j + 1; 
  }
}




/**
 * Generate profile norms
*/
void gen_unorms()
{
  unorm = (fixed_t<2 * param>*)malloc(nprofiles * sizeof(fixed_t<2 * param>));
  _128_unorm = (fixed_t<2 * param, uint128_t>*)malloc(nprofiles * sizeof(fixed_t<2 * param, uint128_t>));
 
  for (size_t j = 0; j < nprofiles; ++j)
  {
    while (norm_progress[norm_step::prenorm_profiles_in] < j + 1 || norm_progress[norm_step::prenorm_profiles_gen] < j + 1 ) 
    { 

 
      std::this_thread::yield(); 
    }
    
    unorm[j] =  (dotprod_recv[j].s +  dot(prenorm_uprofiles[j], prenorm_uprofiles[j] + blinded_prenorm_profile[j]) - dot(blinded_prenorm_profile[j], dotprod_recv[j].u_blind)); 
    
    _128_unorm[j] =  (dotprod_recv[j]._128_s +  dot(_128_prenorm_uprofiles[j], _128_prenorm_uprofiles[j] + _128_blinded_prenorm_profile[j]) 
                      - dot(_128_blinded_prenorm_profile[j], dotprod_recv[j]._128_u_blind)); 

    target[j] =  mod(_128_unorm[j].val >> (3*precision-6), max_target); 

    _128_reduced[j] = (__int128(uint64_t(_128_unorm[j].val))) >> ((3*precision - param));
    
    
    unorm_msb[j] =  (uint64_t(_128_unorm[j].val) & msb) > 0;
    
    if(unorm_msb[j] == 1)
    {
      _128_reduced[j] += ( ~(~0ULL >> (3*precision - param)));
    }

    arc4random_buf(&unorm_msb_this[j], sizeof(uint64_t));

    unorm_msb_written[j] = unorm_msb[j] - unorm_msb_this[j]; 


    norm_progress[norm_step::unorms_out] = j + 1;

  }
}


  /**
   Writes msb of the norm
   @param j
   @param io_context
   @param sout is the socket used to send the unorm msbs
  */
  void write_unorm_msb(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
  {
    while(norm_progress[norm_step::unorms_out] < j + 1)
    {
      //std::cout << "5\n"; 
      std::this_thread::yield();
    }

     async_write(sout, boost::asio::buffer(&unorm_msb_written[j], sizeof(uint64_t)),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      { 
          if(!ec)
          {
             if(j + 1 < nprofiles)
             { 
               write_unorm_msb(j + 1, io_context, sout);
             }
          }          
          else
          {
            write_unorm_msb(j , io_context, sout);
          }
      
      }); 
     
  }  

  /**
   Reads msb of the norm
   @param sin is the socket used to receive the msb of the norm
  */
  void read_unorm_msb(tcp::socket& sin)
  {
    for(size_t j = 0; j < nprofiles; ++j)
    {
      read(sin, boost::asio::buffer(&unorm_msb_recv[j], sizeof(uint64_t)));

      norm_progress[norm_step::msb_in] = j + 1;
    }
  }

  
  /**
   Blinds the msb of the unorm
  */
  void blind_unorm_msb()
  {

    for(size_t j = 0; j < nprofiles; ++j)
    {
        while(norm_progress[norm_step::unorms_out] < j + 1 || norm_progress[norm_step::norm_bvrs_in] < j + 1)
        { 
          std::this_thread::yield();
        }


         blinded_unorm_msb[j] = unorm_msb_this[j] + norm_bvrs_recv[j].blind;

         norm_progress[norm_step::blinded_msb_gen] = j + 1;

    }
    
  }

  
  /**
    * Writes the blinded msb of the norm to Pb
    * @param j
    * @param io_context
    * @param sout is the socket to write
  */
  void write_blinded_unorm_msb(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
  {
    while(norm_progress[norm_step::blinded_msb_gen] < j + 1)
    { 
      std::this_thread::yield();
    }

    

     async_write(sout, boost::asio::buffer(&blinded_unorm_msb[j], sizeof(uint64_t)),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      { 
          if(!ec)
          {
             if(j + 1 < nprofiles){ 
              write_blinded_unorm_msb(j + 1, io_context, sout);
            }
          }          
          else
          {
  
            write_blinded_unorm_msb(j , io_context, sout);
          }
      
      }); 

     
  }

  void read_blinded_unorm_msb(tcp::socket& sin)
  {
    for(size_t j = 0; j < nprofiles; ++j)
    {
      read(sin, boost::asio::buffer(&blinded_unorm_msb_recv[j], sizeof(uint64_t)));
      norm_progress[norm_step::blinded_msb_in] = j + 1;
    }
  }



  /**
   Generates msb product
  */
  void gen_msb_prod()
  {
    for(size_t j = 0; j < nprofiles; ++j)
    {
      while(  norm_progress[norm_step::norm_bvrs_in] < j + 1   ||
              norm_progress[norm_step::msb_in] < j + 1         || 
              norm_progress[norm_step::blinded_msb_in] < j + 1 || 
              norm_progress[norm_step::blinded_msb_gen] < j + 1)
     {
 
      std::this_thread::yield();
     }


     #if(PARTY == 0)  
      msb_prod[j] =  (unorm_msb_recv[j] * unorm_msb_written[j]) + ((unorm_msb_written[j] + unorm_msb_this[j]) * blinded_unorm_msb_recv[j]) - (norm_bvrs_recv[j].blind  * unorm_msb_recv[j]) +   norm_bvrs_recv[j].bb_u;
     #endif

     #if(PARTY == 1)  
      msb_prod[j] =  (unorm_msb_written[j] * blinded_unorm_msb_recv[j])   - (norm_bvrs_recv[j].blind * unorm_msb_recv[j]) - (norm_bvrs_recv[j].blind * blinded_unorm_msb_recv[j]);
     #endif


      norm_progress[norm_step::msb_prod_gen] = j + 1;
    }
  }

  void add_offset(size_t j, boost::asio::io_context& io_context,  tcp::socket& sout)
  {
 
      while(norm_progress[norm_step::msb_prod_gen] < j + 1)
      {
        std::this_thread::yield();
      }
 
 
      async_write(sout, boost::asio::buffer(&msb_prod[j], sizeof(uint64_t)),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      { 
          if(!ec)
          {
             if(j + 1 < nprofiles){ 
              add_offset(j + 1, io_context, sout);
            }
          }          
          else
          {
  
            add_offset(j , io_context, sout);
          }
      
      }); 
 
  }
  
  void read_offset(tcp::socket& sin)
  { 

      for(size_t j = 0; j < nprofiles; ++j)
      {
        while(norm_progress[norm_step::msb_prod_gen] < j + 1) 
        {
          std::this_thread::yield();
        }
     
         read(sin, boost::asio::buffer(&msb_read[j], sizeof(uint64_t)));      
  
         norm_progress[norm_step::offset_in] = j + 1; 
      }
  } 
 

  void subtract_offset()
  {
    for(size_t j = 0; j < nprofiles; ++j)
    {
       
        while(norm_progress[norm_step::offset_in] < j + 1)
        { 
           std::this_thread::yield();
        }
       

        msb_read[j] = msb_read[j] + msb_prod[j];
            
        #if (PARTY == 1)
         reduced[j] -=  msb_read[j] * (~(~0ULL >> (3 * precision - param)));
        #endif
        
        reduced[j] = _128_reduced[j];

        norm_progress[norm_step::offset_subtracted] = j + 1; 

    }
  }


  void read_norm_bvrs(AES_KEY& key, tcp::socket& sin1, tcp::socket& sin2)
  {
    
    __m128i seed;

    read(sin1, boost::asio::buffer(&seed, sizeof(__m128i)));
    
    PRG(key, seed, (__m128i *)bvrs_created, (nprofiles * sizeof(bvrs)) / sizeof(__m128i));

    for(size_t j = 0; j < nprofiles; ++j)
    {
      read(sin2, boost::asio::buffer(&norm_bvrs_recv[j], sizeof(norm_bvrs_recv[j])));

      norm_progress[norm_step::norm_bvrs_in] = j + 1; 
    }
  }
 

  void generate_dpfkeys()
  {  
    for(size_t j = 0; j < nprofiles; ++j)
    {

      while(norm_progress[norm_step::unorms_out] < j + 1 || norm_progress[norm_step::offset_subtracted] < j + 1 ) 
      { 
         std::this_thread::yield();
      }
      
      AES_KEY prgkey;
      
      AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &prgkey);
       
      //std::cout << "target[" << j << "] = " << target[j] << std::endl;

      auto [dpfkey_b, dpfkey_other] = dpf_key<leaf_t, node_t, prgkey_t>::gen(prgkey, max_target, target[j]);
 
      keys_b.push_back(dpfkey_b);
      
      keys_other.push_back(dpfkey_other);

      for(size_t d = 0; d < depth; ++d) dpfkey_written_other[j].cw[d] = keys_other[j].cw[d];
      dpfkey_written_other[j].root = keys_other[j].root;             
      dpfkey_written_other[j].finalizer = keys_other[j].finalizer;


      for(size_t d = 0; d < depth; ++d) dpfkey_written_b[j].cw[d] = keys_b[j].cw[d];
      dpfkey_written_b[j].root = keys_b[j].root; 
      dpfkey_written_b[j].finalizer = keys_b[j].finalizer;

      norm_progress[norm_step::dpf_gen] = j + 1; 
    }
  }
 


 
    void send_share_to_P_other(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
 
       while(norm_progress[norm_step::unorms_out] < j + 1) 
       {  
         std::this_thread::yield(); 
       }
 
        async_write(sout, boost::asio::buffer(&_128_unorm[j].val, sizeof(uint128_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
               if(j + 1 < nprofiles) 
               {
                send_share_to_P_other(j + 1, io_context, sout);
               }
      
            }          
            else
            {               
              send_share_to_P_other(j , io_context, sout);
            }
        
        });
    }


    void send_reduced_to_P_other(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
       
     while(norm_progress[norm_step::offset_subtracted] < j + 1 ) 
     { 
       /// << "11\n"; 
       std::this_thread::yield(); 
     }
    
      async_write(sout, boost::asio::buffer(&reduced[j], sizeof(uint64_t)),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      
      {       
          if(!ec)
          {
            
             if(j + 1 < nprofiles) 
             {
              send_reduced_to_P_other(j + 1, io_context, sout);
             }
    
          }          
          else
          {               
            send_reduced_to_P_other(j , io_context, sout);
          }
      
      });
    }
 
    void compute_s()
    {
      
      for(size_t j = 0; j < nprofiles; ++j)
      {

        while(norm_progress[norm_step::dpfs_in] < j + 1 || norm_progress[norm_step::dpf_gen] < j + 1 || 
              norm_progress[norm_step::unorms_out] < j + 1 || norm_progress[norm_step::offset_subtracted] < j + 1 ) 
        { 
           std::this_thread::yield(); 
        }
      
        query_recv[j] = (uint64_t*)std::aligned_alloc(alignof(__m128i), keys_recv[j].full_bytes());
        
        keys_recv[j].evalfull((bool*)query_recv[j]);

       
        auto [m, b, cnt0] = isqrt_coeffs(query_recv[j], target[j]);
        cnt[j] = cnt0;
       
        free(query_recv[j]); 
      
        #if(PARTY == 1)
          b = 0;
        #endif
 
        long double int_, frac_ = std::modf(m, &int_);
        
 
        
        s[j] = uint64_t(b) + (uint64_t(-int_) * reduced[j]) + (uint64_t(-frac_*scale1_5) * reduced[j]) / scale1_5;
 
        
        norm_progress[norm_step::s_compute] = j + 1;
      }
    }
  

 
    void send_dpfs_to_P_b(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
 
         while(norm_progress[norm_step::dpf_gen] < j + 1) 
          {
           std::this_thread::yield(); 
          }
            
            // std::cout << "-> root written: " << dpfkey_written_b[j].root.mX[0] << " " << dpfkey_written_b[j].root.mX[1] << std::endl; 
            async_write(sout, boost::asio::buffer(&dpfkey_written_b[j], sizeof(dpfkey_written_b[j])),
            [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
            
            {       
                if(!ec)
                {                    
                 if(j + 1 < nprofiles) 
                 {
                  send_dpfs_to_P_b(j + 1, io_context, sout);
                 }
                }          
                else
                {               
                  send_dpfs_to_P_b(j , io_context, sout);
                }
            
            }); 
    } 
 
    void send_dpfs_to_P_other(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
  
           while(norm_progress[norm_step::dpf_gen] < j + 1) 
            {
              std::this_thread::yield(); 
            }
            
            //std::cout << " -> " << dpfkey_written_other[j].root.mX[0] << " " << dpfkey_written_other[j].root.mX[1] << std::endl;

            async_write(sout, boost::asio::buffer(&dpfkey_written_other[j], sizeof(dpfkey_written_other[j])),
            [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
            
            {       
                if(!ec)
                {
                  
                   if(j + 1 < nprofiles) 
                   {
                    send_dpfs_to_P_other(j + 1, io_context, sout);
                   }
          
                }          
                else
                {               
                  send_dpfs_to_P_other(j , io_context, sout);
                }
            
            }); 
 
    } 
 

 
    void read_dpfs_from_P_b(tcp::socket &sin)
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
        size_t nitems_; 
        node_t root_;
        std::vector<node_t> cw_(depth);
        std::array<node_t, nodes_per_leaf> finalizer_;

        AES_KEY prgkey;

        AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &prgkey);

        read(sin, boost::asio::buffer(&dpfkey_read[j], sizeof(dpfkey_read[j])));
        //std::cout << "-> root read: " << dpfkey_read[j].root.mX[0] << " " << dpfkey_read[j].root.mX[1] << std::endl; 
        
        nitems_ = max_target;// 1ULL << 16; //dpfkey_read[j].nitems; 

        root_ = dpfkey_read[j].root;
        
        for(size_t d = 0; d < depth; ++d) cw_[d] = dpfkey_read[j].cw[d];
        
        finalizer_ = dpfkey_read[j].finalizer;
        
        dpf_key<leaf_t, node_t, prgkey_t> dpfkey(std::move(nitems_),  std::move(root_), std::move(cw_), std::move(finalizer_), std::move(prgkey));

        keys_recv.push_back(dpfkey);

        norm_progress[norm_step::dpfs_in] = j + 1;
      }
    }
 
    void recv_s_from_P_other(tcp::socket& sin)
    {   
      for(size_t j = 0; j < nprofiles; ++j)
      {   
       read(sin, boost::asio::buffer(&s_other[j], sizeof(uint64_t)));  
       norm_progress[norm_step::s_other_in] = j + 1; 
      }
    }


     void generate_sshares()
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
          while(norm_progress[norm_step::s_compute] < j + 1 || norm_progress[norm_step::s_other_in] < j + 1) 
          { 
            // std::cout << "15\n"; 
            std::this_thread::yield(); 
          }

         s_s_other[j] = s[j] - s_other[j]; // s - s_other          
         
         arc4random_buf(&sshare0[j], sizeof(uint64_t));
         
         sshare1[j] = s_s_other[j] - sshare0[j];

         norm_progress[norm_step::sshare_gen] = j + 1;
       }
    }

    void compute_z()
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
         while(norm_progress[norm_step::s_compute] < j + 1 || 
               norm_progress[norm_step::norm_bvrs_in] < j + 1 || 
               norm_progress[norm_step::blind_sshare_in] < j + 1 || 
               norm_progress[norm_step::blind_cnt_in] < j + 1 ||
               norm_progress[norm_step::sshare_gen] < j + 1) 
         { 
           // std::cout << "14\n"; 
             std::this_thread::yield(); 
         }

        z[j] = sshare0[j] * (cnt[j] + blinded_cnt_recv[j]) - (bvrs_created[j].s_blind0 * blinded_sshare_recv[j]) + norm_bvrs_recv[j].s_gamma;

        norm_progress[norm_step::z_compute] = j + 1;
      }
    }
 



    void write_sshare(size_t j, boost::asio::io_context& io_context, tcp::socket &sout)
    {
        while(norm_progress[norm_step::sshare_gen] < j + 1) 
          { 
             std::this_thread::yield(); 
          }
   
          
        async_write(sout, boost::asio::buffer(&sshare1[j], sizeof(uint64_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
               if(j + 1 < nprofiles) 
               {
                write_sshare(j + 1, io_context, sout);
               }
      
            }          
            else
            {               
              write_sshare(j , io_context, sout);
            }
        
        });
         
       
    } 
    

    void gen_blind_sshares()
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
        while(norm_progress[norm_step::norm_bvrs_in] < j + 1 || norm_progress[norm_step::sshare_gen] < j + 1) 
        { 
          std::this_thread::yield(); 
        }

        blind_sshare[j]  = sshare0[j] + bvrs_created[j].s_blind1;  

        norm_progress[norm_step::blind_sshare_gen] = j + 1;
      }
    }


    void write_blind_sshare(size_t j,  boost::asio::io_context& io_context, tcp::socket &sout)
    {
   
        while(norm_progress[norm_step::blind_sshare_gen] < j + 1) 
        { 
          std::this_thread::yield(); 
        }

        async_write(sout, boost::asio::buffer(&blind_sshare[j], sizeof(uint64_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
               if(j + 1 < nprofiles) 
               {
                write_blind_sshare(j + 1, io_context, sout);
               }
      
            }          
            else
            {               
              write_blind_sshare(j , io_context, sout);
            }
        
        });
       
    }
    


    void generate_blind_cnt()
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
        while(norm_progress[norm_step::s_compute] < j + 1 || norm_progress[norm_step::norm_bvrs_in] < j + 1) 
        { 
          std::this_thread::yield(); 
        }
   
        blind_cnt[j]  = cnt[j] + bvrs_created[j].s_blind0; 

        norm_progress[norm_step::blind_cnt_gen] = j + 1;
      }
    }

    void write_blind_cnt(size_t j, boost::asio::io_context& io_context, tcp::socket &sout)
    {
 
        while(norm_progress[norm_step::blind_cnt_gen] < j + 1) 
        { 
          std::this_thread::yield(); 
        }       

        async_write(sout, boost::asio::buffer(&blind_cnt[j], sizeof(uint64_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        {       
            if(!ec)
            {
              
               if(j + 1 < nprofiles) 
               {
                write_blind_cnt(j + 1, io_context, sout);
               }
      
            }          
            else
            {               
              write_blind_cnt(j , io_context, sout);
            }
        
        }); 
    }
    


    void read_blind_sshare(tcp::socket &sin)
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
        read(sin, boost::asio::buffer(&blinded_sshare_recv[j], sizeof(uint64_t)));
        norm_progress[norm_step::blind_sshare_in] = j + 1;
      }
    }


    void read_blind_cnt(tcp::socket &sin)
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
       read(sin, boost::asio::buffer(&blinded_cnt_recv[j], sizeof(uint64_t)));
       norm_progress[norm_step::blind_cnt_in] = j + 1;
      }
    }

  


    void read_z_from_P_other(tcp::socket& sin)
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
          read(sin,  boost::asio::buffer(&zz[j], sizeof(int64_t)));
          norm_progress[norm_step::zz_in] = j + 1;
      }

    }


    void gen_z()
    {
      for(size_t j = 0; j < nprofiles; ++j)
      {
         while(norm_progress[norm_step::z_compute] < j + 1 || norm_progress[norm_step::zz_in] < j +1) 
         {    
          std::this_thread::yield(); 
         } 
       
        z_added[j] =  z[j] + zz[j];

        norm_progress[norm_step::z_gen] = j + 1;
      }
    }

    void write_z_added(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
    {
   
       while(norm_progress[norm_step::z_gen] < j + 1) 
       {  
         std::this_thread::yield(); 
       } 
       
       async_write(sout, boost::asio::buffer(&z_added[j], sizeof(uint64_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
             if(j + 1 < nprofiles) 
             {
              write_z_added(j + 1, io_context, sout);
             }
      
            }          
            else
            {               
              write_z_added(j , io_context, sout);
            }
        
        }); 

        norm_progress[norm_step::z_added_done] = j + 1;
 
   }

   void read_z_added(tcp::socket & sin)
   {
     for(size_t j = 0; j < nprofiles; ++j)
     {

       while(norm_progress[norm_step::z_compute] < j + 1 || norm_progress[norm_step::z_added_done] < j + 1) 
       {
        std::this_thread::yield(); 
       } 
           
      read(sin, boost::asio::buffer(&z_recv[j], sizeof(uint64_t)));
   

      fixed_t<(3 * precision) - param> isqrt =  double(z_recv[j] + z_added[j])/scale1_5;      
 
      for(size_t dim = 0; dim < DIM; ++dim)
      {
        normalized_profile_this[j].f64[dim].val = prenorm_uprofiles[j].f64[dim].val * isqrt.val;
      }
 
      std::cout << "isqrt[" << j << "] = " << isqrt << " = "  << double(z_recv[j] + z_added[j])/scale1_5 << std::endl;
 
      norm_progress[norm_step::norm_done] = j + 1;
    }
   }

   void write_normalized_shares(size_t j, boost::asio::io_context& io_context, tcp::socket & sout)
   { 

      while(norm_progress[norm_step::norm_done] < j + 1)
      {
        std::this_thread::yield();
      }
 
       async_write(sout, boost::asio::buffer(&normalized_profile_this[j], sizeof(profile<3 * precision>)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
               if(j + 1 < nprofiles) 
               {
                write_normalized_shares(j + 1, io_context, sout);
               }
      
            }          
            else
            {               
              write_normalized_shares(j , io_context, sout);
            }
        
        });

      norm_progress[norm_step::normalized_out] = j + 1;
   }

   void read_normalized_shares(tcp::socket & sin)
   {
      for(size_t j = 0; j < nprofiles; ++j)
      {      
          while(norm_progress[norm_step::norm_done] < j + 1)
          {
            std::this_thread::yield();
          }
          read(sin,  boost::asio::buffer(&normalized_profile_other[j], sizeof(profile<3 * precision>)));

          std::cerr << "reconstructed " << j << " : -> " << (normalized_profile_this[j] + normalized_profile_other[j]) << std::endl;

          norm_progress[norm_step::normalized_in] = j + 1;
      }
   }

int main(int argc, char * argv[])
{  

  keys_recv.reserve(nprofiles);
  keys_b.reserve(nprofiles);
  keys_other.reserve(nprofiles);

  const std::string host1 = (argc < 2) ? "127.0.0.1" : argv[1];
  const std::string host2 = (argc < 3) ? "127.0.0.1" : argv[2];
  const std::string host3 = (argc < 4) ? "127.0.0.1" : argv[3];
  
  std::string ufile = std::string("./updated_uprofiles.") + std::to_string(PARTY);
  // memory map the user profile shares
  int ufd = open(ufile.c_str(), O_RDWR);
  prenorm_uprofiles_ = (profile<3 * (precision/2)> *) mmap(NULL,
  nprofiles * sizeof(profile<3 * (precision/2)>), PROT_READ, MAP_PRIVATE, ufd, 0);

  for (size_t i = 0; i < nprofiles; ++i)
  {
    //prenorm_uprofiles[i] = prenorm_uprofiles_3P[i];
    prenorm_uprofiles[i] = prenorm_uprofiles_[i]; 
    //std::cout << "prenorm_uprofiles[" << i << "] = " << prenorm_uprofiles[i] << std::endl;
    for(size_t d = 0; d < DIM; ++d) _128_prenorm_uprofiles[i].u128[d] = prenorm_uprofiles[i].u64[d];
  }
 

 AES_KEY prgkey;

 AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &prgkey);
  
 
  #include "connections.h"

 
 


  std::thread poller(keep_polling, std::ref(io_context));

  std::thread dotprod_bvrs_reader(read_dotprod_bvrs, std::ref(s2_dotprod));

  std::thread norm_bvrs_reader(read_norm_bvrs, std::ref(prgkey), std::ref(sother_dpf), std::ref(sother_dpf_b));

 
  

  
  std::thread prenorm_profile_blinder(blind_prenorm_profiles);
  
  std::thread profile_reader(read_prenorm_profiles, std::ref(io_context), std::ref(sb_b));
  
  std::thread profile_writer(write_prenorm_profiles,0, std::ref(io_context), std::ref(sb_b));
  
  std::thread unorm_generator(gen_unorms);


  
  std::thread msb_writer(write_unorm_msb, 0, std::ref(io_context), std::ref(sb_c));
  
  std::thread msb_reader(read_unorm_msb, std::ref(sb_c));

  std::thread unorm_msb_blinder(blind_unorm_msb);

  std::thread blinded_unorm_msb_writer(write_blinded_unorm_msb, 0, std::ref(io_context), std::ref(sb_d));
  
  std::thread blinded_unorm_msb_reader(read_blinded_unorm_msb, std::ref(sb_d));
 
  std::thread msb_prod_generator(gen_msb_prod);
 
  std::thread offset_adder(add_offset, 0, std::ref(io_context), std::ref(sb_e));

  std::thread offset_reader(read_offset, std::ref(sb_e));

  std::thread offset_subtractor(subtract_offset);



  std::thread dpf_generator(generate_dpfkeys);






  std::thread share_to_P_other_sender(send_share_to_P_other, 0, std::ref(io_context), std::ref(sother));

  std::thread reduced_to_P_other_sender(send_reduced_to_P_other, 0, std::ref(io_context), std::ref(sother_i));

  std::thread dpfs_to_P_b_sender(send_dpfs_to_P_b, 0, std::ref(io_context), std::ref(sb_f));

  std::thread dpfs_to_P_other_sender(send_dpfs_to_P_other, 0, std::ref(io_context), std::ref(sother_dpf_a));

  std::thread dpf_reader(read_dpfs_from_P_b, std::ref(sb_f));

  std::thread s_computer(compute_s);

  std::thread s_from_P_other_receiver(recv_s_from_P_other, std::ref(sother_a));
  
  std::thread sshares_generator(generate_sshares);

  std::thread sshare_writer(write_sshare, 0,  std::ref(io_context), std::ref(sother_b));

  std::thread blind_sshare_generator(gen_blind_sshares);

  std::thread blind_sshare_sender(write_blind_sshare, 0, std::ref(io_context),  std::ref(sother_c));


  std::thread blind_cnt_generator(generate_blind_cnt);

  std::thread blind_cnt_sender(write_blind_cnt, 0, std::ref(io_context),  std::ref(sother_d));

  std::thread blind_sshare_reader(read_blind_sshare, std::ref(sother_e));


  std::thread blind_cnt_reader(read_blind_cnt, std::ref(sother_f));
  
  std::thread z_computer(compute_z);
 
  std::thread z_from_P_other_reader(read_z_from_P_other, std::ref(sother_g)); 
   
  std::thread z_generator(gen_z);
    
  std::thread z_added_writer(write_z_added, 0, std::ref(io_context), std::ref(sb_g));
  
  std::thread z_added_reader(read_z_added, std::ref(sb_g));

  std::thread normalized_profile_writer(write_normalized_shares, 0, std::ref(io_context), std::ref(sb_h));
  
  std::thread normalized_profile_reader(read_normalized_shares, std::ref(sb_h));
  
  poller.join();
  
  dotprod_bvrs_reader.join();
  
  norm_bvrs_reader.join();

 


  prenorm_profile_blinder.join();

  profile_reader.join();

  profile_writer.join();
  
  unorm_generator.join();
  
  msb_writer.join();
  
  msb_reader.join();

  unorm_msb_blinder.join();
  
  blinded_unorm_msb_writer.join();
  
  blinded_unorm_msb_reader.join();
  
  msb_prod_generator.join();

  offset_adder.join();

  offset_reader.join();

  offset_subtractor.join();

  dpf_generator.join();

  share_to_P_other_sender.join();

  reduced_to_P_other_sender.join();

  dpfs_to_P_b_sender.join();

  dpfs_to_P_other_sender.join();

  dpf_reader.join();

  s_computer.join();

  s_from_P_other_receiver.join();

  sshares_generator.join();

  sshare_writer.join();

  blind_sshare_generator.join();

  blind_sshare_sender.join();

  blind_cnt_generator.join();

  blind_cnt_sender.join();

  blind_sshare_reader.join();

  blind_cnt_reader.join();

  z_computer.join();

  z_from_P_other_reader.join();

  z_generator.join();

  z_added_writer.join();

  z_added_reader.join();

  normalized_profile_writer.join();
  
  normalized_profile_reader.join();

  return 0;
}

 