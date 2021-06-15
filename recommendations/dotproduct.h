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

 uint64_t uprofile[DIM];  // The user profile share
 uint64_t profiles[nrecords + 1][DIM]; // The item profile shares
 uint64_t blinded_user_profiles[DIM];
 uint64_t blinded_item_profiles[nrecords + 1][DIM];
 uint64_t blinded_item_profiles_recv[nrecords + 1][DIM];
 uint64_t blinded_user_profiles_recv[DIM];
 uint64_t profile_blinds[2][nrecords][DIM];
 uint64_t profile_alphas[nrecords];
 uint64_t profile_dots[nrecords];
 block<__m256i>  __profile_dots[nrecords/4];

/**
 * computes the dot-product
*/
 uint64_t dot_prod(uint64_t a[DIM], uint64_t b[DIM])
 {

   uint64_t dot_prod_ = 0;
   for(size_t dim = 0; dim < DIM; ++dim)
   {
      dot_prod_ += (a[dim] * b[dim]);
   }

   return dot_prod_;
 }

/**
 * Reads the cancellation terms from P2 to compute the Du-Attalah dot-product  
 * @param sin is the socket from which P_{b} receives the cancellation terms to compute the dot-product 
*/
void read_profile_alphas_from_P2(tcp::socket &sin)
{
  for(size_t j = 0; j < nrecords; ++j)
  {
     read(sin, boost::asio::buffer(&profile_alphas[j], sizeof(profile_alphas[j])));
     comparison_progress[comparison_step::profile_alphas_in] = j + 1;
  }
}  




 /**
   * Blinds the user profile \u
   * Blinds the item profiles \v[0], \ldots, \v[k] (the profiles of items downloaded by u)
 */
 void blind_profiles()
 {

    for(size_t dim = 0; dim < DIM; ++dim) 
    {
      blinded_user_profiles[dim] = uprofile[dim] + profile_blinds[1][0][dim];
    }

    for(size_t j = 0; j < nrecords; ++j)
    {
       for(size_t dim = 0; dim < DIM; ++dim) 
       {  
         blinded_item_profiles[j][dim] = profiles[j][dim] + profile_blinds[0][j][dim];
       }

      comparison_progress[comparison_step::blinded_profiles_done] = j + 1;
    }
 }


/**
  * P_b writes the blinded user profile to P_{b-1}.
*/
void  write_blinded_user_profiles(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{


         while
              (
               comparison_progress[comparison_step::blinded_profiles_done] < 1
              )
       {
         std::this_thread::yield();
       }
     
     async_write(sout, boost::asio::buffer(&blinded_user_profiles[j],  sizeof(blinded_user_profiles[j])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t)  
      {       
          if(!ec)
          {
             if(j + 1 < DIM) 
             {             
              write_blinded_user_profiles(j + 1, io_context, sout);
             }
          }          
          else
          {               
            write_blinded_user_profiles(j , io_context, sout);
          }     
      }); 
}


/**
 * P_b reads the blinded profile from P_{b-1}
 * @param takes in the incoming socket
*/

void read_blinded_user_profiles(tcp::socket& sin)
{
   for(size_t j = 0; j < DIM; ++j)
   {
     read(sin, boost::asio::buffer(&blinded_user_profiles_recv[j],  sizeof(blinded_user_profiles_recv[j])));
     comparison_progress[comparison_step::blinded_user_profiles_in] = j + 1;
   }
}



/**
 P_b writes the blinded_item_profiles to P_{b-1} 
 @param j
 @param io_context
 @param sout is the socket through blinded item profiles are written.

*/
void  write_blinded_item_profiles(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{
         while(
               comparison_progress[comparison_step::blinded_profiles_done] < j + 1
             )
       {
        std::this_thread::yield();
       }
     
     async_write(sout, boost::asio::buffer(blinded_item_profiles[j], DIM * sizeof(blinded_item_profiles[j][0])),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t)  
      {       
          if(!ec)
          {

             if(j + 1 < nrecords) 
             {             
              write_blinded_item_profiles(j + 1, io_context, sout);
             }
          }          
          else
          {               
            write_blinded_item_profiles(j , io_context, sout);
          }

      
      }); 
}


/**
  * P_b reads the blinded item profiles from P_{b-1}
  * @param sin is the socket with which is used to receive the blinded item profiles

*/

void read_blinded_item_profiles(tcp::socket& sin)
{
   for(size_t j = 0; j < nrecords; ++j)
   {
  
     read(sin, boost::asio::buffer(blinded_item_profiles_recv[j], DIM * sizeof(blinded_item_profiles_recv[j][0])));
      
     comparison_progress[comparison_step::blinded_item_profiles_in] = j + 1;
   }
}


/**
  * P0 and P1 compute the shares of the dot-product between \u and \v[i], where i \in \{0, \ldots, nrecoreds}
*/
void compute_dot_products()
{

  for(size_t j = 0; j < nrecords; ++j)
  {
    while (
           comparison_progress[comparison_step::blinded_item_profiles_in] < j + 1
          )
    {
       std::this_thread::yield();
    }
    profile_dots[j] =   dot_prod(uprofile, blinded_item_profiles_recv[j]) - dot_prod(profile_blinds[0][j], blinded_user_profiles_recv)  + profile_alphas[j];
    comparison_progress[comparison_step::profile_dots_done] = j + 1;
  }
  
   start = std::chrono::steady_clock::now();

}
