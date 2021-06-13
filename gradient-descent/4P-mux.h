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

#include <memory>

void read_mux_bvr(AES_KEY &key, tcp::socket & sin_seed, tcp::socket & sin)
{ 
  
   read(sin_seed, boost::asio::buffer(&seed3, sizeof(__m128i)));
   progress[step::p3_seed_in] = 1;
 
   #if(PARTY == 0)

     std::cout << "seed: " << seed3[0] << " " << seed3[1] << std::endl;
    
    
     const size_t buflen = (nitems) * sizeof(profile<precision>);
    
     PRG(key, seed3, (__m128i *)A, buflen / sizeof(__m128i));
    
     progress[step::mux_seed] = 1;
    
   #endif
 
 
   p3_bvrs = (p3_beavers *) std::aligned_alloc(sizeof(__m256i), nqueries * sizeof(p3_beavers)); 
 
   for(size_t i = 0; i < nqueries; ++i)
   {

     #ifdef VERBOSE
      std::cout << " (p3 beavers) i = " << i << std::endl;
     #endif  
     
     boost::asio::read(sin, boost::asio::buffer(&p3_bvrs[i], p3_beavers::size));
     
     progress[step::bvr3_in] = i + 1;
   }

}



/**
 * Writes blinded items for computing the item-profile shares
*/
void write_blinded_items(const size_t i, boost::asio::io_context& io_context, tcp::socket& sout)
{
      while(progress[step::itm_blind_gen] < i + 1)
      { 
        std::this_thread::yield();
      }
      
     async_write(sout, boost::asio::buffer(&A[i], sizeof(profile<precision>)),
            [i, &io_context, &sout](boost::system::error_code ec, std::size_t) 
            
            {                    
                if(!ec)
                {
                   if(i + 1 < nitems)
                   { 
                    
                     #ifdef VERBOSE
                      std::cout << "write_blinded_items[" << i << "] < " << progress[step::itm_blind_gen] << std::endl;
                     #endif
                     
                     write_blinded_items(i + 1, io_context, sout);
                   }
                }          
                else
                {
                  write_blinded_items(i , io_context, sout);
                }
            
            }); 

     progress[step::itm_blind_out] = i + 1;

}


 
/**
 * Blinds the item profiles
*/
void generate_blinded_items()
{

   for(size_t i = 0; i < nitems; ++i)
   {
      while (progress[step::mux_seed] < 1) 
      { 
         #ifdef VERBOSE    
          std::cout << " stuck at write_blinded_items: " << i << std::endl;
         #endif
          //std::cout << "21\n";  
         std::this_thread::yield(); 
      }
 
      
    A[i] = iprofiles[i] + A[i];
     
    progress[step::itm_blind_gen] = i + 1;
   
   } 
}


/**
 * Reads blinded items
*/
void read_blinded_items(tcp::socket & sin)
{
 
   itmBB = (profile<precision> *) std::aligned_alloc(sizeof(__m256i), nitems * sizeof(profile<precision>));
    
   for (size_t i = 0; i < nitems; ++i)
   { 
      while (progress[step::mux_seed] < 1) 
      {
 
         std::this_thread::yield(); 
 
      }
    
      read(sin, boost::asio::buffer(&itmBB[i], sizeof(profile<precision>)));
       
      itmBB[i] = iprofiles[i] + itmBB[i];
    
      progress[step::itm_blind_in] = i + 1;
   }
  
}



void write_mux_norm(size_t i, boost::asio::io_context & io_context, tcp::socket & sout)
{
   

  pm1 = (fixed_t<2*precision>*)std::aligned_alloc(sizeof(__m256i),nqueries * sizeof(fixed_t<2*precision>));
  
  itm = (profile<precision> *) std::aligned_alloc(sizeof(__m256i), nqueries * sizeof(profile<precision>));

 

  for(size_t i = 0; i < nqueries; ++i)
  {

   while(
          progress[step::bvr2_sanity_in] < i + 1 ||
          progress[step::bvr3_in] < i + 1        ||
          progress[step::itm_blind_in] < nitems  ||
          progress[step::itm_blind_gen] < nitems || 
          progress[step::evalfull_done_]  < i + 1
          ) 

        {
 
          std::this_thread::yield(); 
        }

    //if(sanity_vector_demux[i] !=  2)
    {
     dpf_mux2<nitems,PARTY>(aeskey, queries[i].dpfkey, itmBB, itm[i], S[i], T[i], pm1[i], i);
    }
    // else
    // {
    //   std::cerr << i << ": " << "dpf_mux" << std::endl;
    // }
  
  
    #if(PARTY == 1)
    itm[i] = itm[i] +  p3_bvrs[i].A_m;  
    #endif

    #if(PARTY == 0)
    itm[i] = itm[i] + bvr2_sanity[i].A_m;// + blind_s;   
    #endif
    
    progress[step::itm_mux_out] = i + 1;
  }
 
  free(itmBB);
  free(p3_bvrs);
  #if(PARTY == 1) 
    free(A);
  #endif
}


 