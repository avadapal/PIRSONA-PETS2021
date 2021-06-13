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



/**
 * blinds the I_profile msbs
*/
void blind_I_profile_msbs(AES_KEY& prgkey)
{
  while(progress[step::precision_seed_in] < 1) 
  {
    std::this_thread::yield();
  }
 
  PRG(prgkey, precision_seed, (__m128i *)I_profile_msb_blind, (sizeof(uint64_t) * total_bucket_elements * DIM) / sizeof(__m128i));
 
  for(size_t j = 0; j < total_bucket_elements; ++j)
  {
     while(progress[step::i_msb_in] < j + 1) 
     {
       std::this_thread::yield();
     }
   
    for(size_t dim = 0; dim < DIM; ++dim)
    { 
      I_profile_msbs_blinded[j][dim] = I_profile_msbs_this[j][dim] + I_profile_msb_blind[j][dim];
    }

    progress[step::i_profile_msbs_blinded] = j + 1;
  }
}


/**
 * Writes the blinded profile msbs to Pb
*/
void write_blinded_I_profile_msbs(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{

    while(progress[step::i_profile_msbs_blinded] < j + 1)
    {
      std::this_thread::yield();
    }

     async_write(sout, boost::asio::buffer(I_profile_msbs_blinded[j], DIM * sizeof(uint64_t)),
          [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
          
          {       
              if(!ec)
              {
                
                 if(j + 1 < total_bucket_elements)
                 {
                   write_blinded_I_profile_msbs( j + 1, io_context, sout);
                 }
              }          
              else
              {
                 write_blinded_I_profile_msbs(j , io_context, sout);
              }
          
          });  
}



/**
 * Reads the blinded profile msbs from Pb
*/
void read_blinded_I_profile_msbs(tcp::socket& sin)
{
  for(size_t j = 0; j < total_bucket_elements; ++j)
  {
    read(sin, boost::asio::buffer(I_profile_msbs_blinded_read[j], DIM * sizeof(uint64_t)));

    progress[step::blinded_i_msbs_in] = j + 1;
  }
}


void gen_msb_prod()
{
  for(size_t j = 0; j < total_bucket_elements; ++j)
  {
    while(progress[step::precision_seed_in] < 1 || progress[step::i_profile_msbs_blinded] < j + 1 || progress[step::blinded_i_msbs_in] < j + 1
          || progress[step::precision_beaver_in] < j + 1 || progress[step::iprofile_msbs_in] < j + 1)
    {
      std::this_thread::yield();
    }

    for(size_t dim = 0; dim < DIM; ++dim)
    {

      #if(PARTY == 0)  
         msb_prod[j][dim] =  (I_profile_msbs_read[j][dim] * I_profile_msbs_written[j][dim]) + ((I_profile_msbs_written[j][dim] + I_profile_msbs_this[j][dim]) * I_profile_msbs_blinded_read[j][dim]) 
                                  - (I_profile_msb_blind[j][dim] * I_profile_msbs_read[j][dim]) + I_profile_msb_alpha[j][dim];
      #endif

      #if(PARTY == 1)  
       msb_prod[j][dim] =  (I_profile_msbs_written[j][dim] * I_profile_msbs_blinded_read[j][dim])  - (I_profile_msb_blind[j][dim] * I_profile_msbs_read[j][dim])
                            - (I_profile_msb_blind[j][dim] * I_profile_msbs_blinded_read[j][dim]);
      #endif

      assert( msb_prod[j][dim] != 0);
    
    }

    progress[step::msb_prod_gen] = j + 1;
  }
}


void write_msb_prod(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{

  while(progress[step::msb_prod_gen] < j + 1)
  {
    std::this_thread::yield();
  }



 
   async_write(sout, boost::asio::buffer(msb_prod[j], DIM * sizeof(uint64_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {
              
               if(j + 1 < total_bucket_elements)
               {
                 write_msb_prod( j + 1, io_context, sout);
               }
            }          
            else
            {
               write_msb_prod(j , io_context, sout);
            }
        
        });  
}

void read_msb_prod(tcp::socket& sin)
{
  for(size_t j = 0; j < total_bucket_elements; ++j)
  {
     while(progress[step::msb_prod_gen] < j + 1)
     {
       std::this_thread::yield();
     }

    read(sin, boost::asio::buffer(&msb_prod_recv[j], DIM * sizeof(uint64_t)));

    for(size_t dim = 0; dim < DIM; ++dim) 
     {
        #ifdef VERBOSE
         std::cout << " msbsum-> " << (msb_prod_recv[j][dim] + msb_prod[j][dim]) << " = " << msb_prod_recv[j][dim] <<  " + " << msb_prod[j][dim] << std::endl;
        #endif
         
        if(msb_prod_recv[j][dim] + msb_prod[j][dim] != 0)
        {
          std::cerr << " msbsum-> " << (msb_prod_recv[j][dim] + msb_prod[j][dim]) << " = " << msb_prod_recv[j][dim] <<  " + " << msb_prod[j][dim] << std::endl;
        }
     }
  }
}


/**
 * Reads beavers for precision reduction
 * @param sin is the socket to read from P2
*/
void read_precision_reduction_bvrs(tcp::socket& sin)
{

  read(sin, boost::asio::buffer(&precision_seed, sizeof(__m128i))); 
  std::cerr << "precision_seed = " << precision_seed[0] << " " << precision_seed[1] << std::endl;
  progress[step::precision_seed_in] = 1; 
 
  for(size_t j = 0; j < total_bucket_elements; ++j)
  {
     read(sin, boost::asio::buffer(&I_profile_msb_alpha[j], DIM * sizeof(uint64_t)));

     progress[step::precision_beaver_in] = j + 1;
  }
}



/**
 * writes the msbs of the profiles
*/
void write_I_profile_msbs(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{
     while(progress[step::i_msb_in] < j + 1)
     {
        std::this_thread::yield();
     }

     if(j % bucket_size_ == 0 && precision == 16)
     {
       io_context.poll();
     } 

     async_write(sout, boost::asio::buffer(I_profile_msbs_written[j], DIM * sizeof(uint64_t)),
          [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
          
          {       
              if(!ec)
              {
                
                 if(j + 1 < total_bucket_elements)
                 {
                   write_I_profile_msbs( j + 1, io_context, sout);
                 }
              }          
              else
              {
                 write_I_profile_msbs(j , io_context, sout);
              }
          
          });  
}


/**
 * reads msbs for precision reduction
 * @param sin is the socket to read msbs 
*/
void read_I_profile_msbs(tcp::socket& sin)
{

  std::cerr << "reading\n";

  size_t items_read = 0;
  for(size_t j = 0; j < total_bucket_elements; ++j)
  { 
     while(progress[step::i_msb_in] < j + 1)
     {
        std::this_thread::yield();
     }

    read(sin, boost::asio::buffer(I_profile_msbs_read[j], DIM * sizeof(uint64_t)));

    progress[step::iprofile_msbs_in] = j + 1; 
  }
  

} 


/**
 * Computes the final precision reduced vectors
*/
void adjust_offset()
{
  for(size_t j = 0; j < total_bucket_elements; ++j)
  {
    while(
      progress[step::iprofile_msbs_in] < j + 1 
          || progress[step::msb_prod_gen] < j + 1
          )
    {
      std::this_thread::yield();
    }

   for(size_t dim = 0; dim < DIM; ++dim)
    {
     
      iprofiles_reduced[(j % nitems)].u64[dim] -= msb_prod[j][dim] * ( ~(~0ULL >> (3*precision - (3 * precision/2))) );
     
      if(I_profile_msbs_read[j][dim] == 1 && I_profile_msbs[j][dim] == 1)
      {
        std::cerr << "Both One!\n\n " << (j % nitems) << " " << (j / bucket_size_) << std::endl;  
      }
    }
    
    progress[step::offset_iprofiles] = j + 1; 
  }
}