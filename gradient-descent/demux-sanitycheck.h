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


void demux_sanity_check(AES_KEY& prgkey)
{
  
  dp_sanity = (uint64_t *) std::aligned_alloc(sizeof(__m256i), nqueries * sizeof(uint64_t));  

  profile<3*precision>  tCW; 

  profile<3*precision>  zero_profile;
  
  for(size_t dim = 0; dim < DIM; ++dim) zero_profile.u64[dim] = 0;

  for(size_t j = 0; j < nqueries; ++j)
  {
      while (progress[step::evalfull_done_] < j + 1 || progress[step::leaf_expanded_] < j + 1 || progress[step::bvr3_in] < j + 1 || progress[bvr2_sanity_in] < j + 1) 
      { 
        std::this_thread::yield(); 
      }
       
      uint64_t  RRtCW = 0.0;
      
      for(size_t k = 0; k < nitems; ++k)
      { 
       
        if(T[j][k])  tCW  = queries[j].dpfkey.leaf;
        
        if(!T[j][k]) tCW  =  zero_profile; 
        
        tCW = expanded_leaf[j][k] - tCW - tCW;  
        
        for(size_t dim = 0; dim < DIM; ++dim) RRtCW += tCW.u64[dim];
      }

      #if (PARTY == 0)
        dp_sanity[j] = RRtCW + bvr2_sanity[j].RRtCW;// - p3_bvrs[j].sanity_blinds;    
      #endif

      #if (PARTY == 1)
        dp_sanity[j] = RRtCW + p3_bvrs[j].RRtCW;// - bvr2[j].sanity_blinds;  
      #endif

      progress[step::demux_sanity_gen] = j + 1;
  }
}


 void write_demux_sanity_dp(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
 {
    while(progress[step::demux_sanity_gen] < j + 1)
    {
      std::this_thread::yield();
    }

     async_write(sout, boost::asio::buffer(&dp_sanity[j], sizeof(uint64_t)),
            [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
            
            {       
                if(!ec)
                {                   
                   if(j + 1 < nqueries) 
                   {
                    write_demux_sanity_dp(j + 1, io_context, sout);
                   }
                }          
                else
                {                
                  write_demux_sanity_dp(j , io_context, sout);
                }
            
            }); 

 }

 void read_demux_sanity_dp(tcp::socket& sin)
 {
   dp_sanity_recv = (uint64_t *) std::aligned_alloc(sizeof(__m256i), nqueries * sizeof(uint64_t ));
   //uint64_t dp_sanity_recv;
   size_t demux_verified = 0;
   for(size_t j = 0; j < nqueries; ++j)
   {
      while(progress[step::demux_sanity_gen] < j + 1)
      {
       std::this_thread::yield();
      }

      boost::asio::read(sin, boost::asio::buffer(&dp_sanity_recv[j], sizeof(uint64_t)));

      if(uint64_t(dp_sanity_recv[j]) == uint64_t(dp_sanity[j])) 
      {
        sanity_vector_demux[j] = 1;
        ++demux_verified;
      }
 
      
      progress[step::demux_sanity_done] = j + 1;
   }

   std::cerr << "demux_verified = " << demux_verified << std::endl;
 }