   size_t s[nprofiles]; 
   int64_t cnt[nprofiles] = {0};
   uint64_t blind_cnt[nprofiles]; 
   uint64_t sshare_recv[nprofiles];
   uint64_t blind_sshare[nprofiles];
   uint64_t blinded_cnt_recv[nprofiles];
   uint64_t blinded_sshare_recv[nprofiles];
   int64_t z[nprofiles];
   uint128_t X[nprofiles]; 
   uint64_t reduced[nprofiles];
 
   uint32_t target[nprofiles];






void send_blinds_to_P_b(size_t j, boost::asio::io_context &io_context, tcp::socket& sout)
{


		
      while(norm_progress[norm_step::blinds_created] < j + 1) 
      { 
        std::this_thread::yield(); 
      }

     boost::asio::async_write(sout, boost::asio::buffer(&norm_bvrs_b[j], sizeof(norm_bvrs_b[j])),
           [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      { 
          if(!ec)
          {
           if (j+1 < nprofiles) 
            {   
              send_blinds_to_P_b(j+1, io_context, sout); 
            }
          }
          else
          {
            send_blinds_to_P_b(j, io_context, sout);
          } 
      });

    norm_progress[norm_step::norm_bvrs_PB_out] = j + 1;
 }


  void send_blinds_to_P_other(size_t j, boost::asio::io_context &io_context, tcp::socket& sout)
  {	 
       
     while(norm_progress[norm_step::blinds_created] < j + 1) 
     { 
       std::this_thread::yield(); 
     }
        	
     boost::asio::async_write(sout, boost::asio::buffer(&norm_bvrs_other[j], sizeof(norm_bvrs_other[j])),
          [j, &io_context, &sout](boost::system::error_code ec, std::size_t)
          { 
            if(!ec)
            {
             if (j+1 < nprofiles) 
               {   
                send_blinds_to_P_other(j+1, io_context, sout); 
               }
            }
            else
            {
              send_blinds_to_P_other(j, io_context,  sout);
            } 

        });

    norm_progress[norm_step::norm_bvrs_P_other_out] = j + 1;
  }

void create_blinds(AES_KEY & key,  tcp::socket& s_b, tcp::socket& s_other)
{

   
  __m128i seed[2];

   arc4random_buf(seed, sizeof(__m128i) * 2);   
   write(s_b, boost::asio::buffer(&seed[0], sizeof(__m128i))); 
   write(s_other, boost::asio::buffer(&seed[1], sizeof(__m128i))); 

   PRG(key, seed[0], (__m128i *)bvrs_b, (nprofiles * sizeof(bvrs)) / sizeof(__m128i));
   PRG(key, seed[1], (__m128i *)bvrs_other, (nprofiles * sizeof(bvrs)) / sizeof(__m128i));

    for(size_t j = 0; j < nprofiles; ++j)
    {
          uint64_t alpha;
          arc4random_buf(&alpha, sizeof(uint64_t));
          //std::cout << "alpha = " << alpha << std::endl;
          norm_bvrs_b[j].r_gamma     = bvrs_b[j].r_blind0 * bvrs_other[j].r_blind1;// + alpha;
          norm_bvrs_other[j].r_gamma = bvrs_b[j].r_blind1 * bvrs_other[j].r_blind0;// - alpha;

          arc4random_buf(&alpha, sizeof(uint64_t));
         // std::cout << "alpha = " << alpha << std::endl;
          norm_bvrs_b[j].s_gamma = bvrs_b[j].s_blind0 * bvrs_other[j].s_blind1;// + alpha;
          norm_bvrs_other[j].s_gamma = bvrs_b[j].s_blind1 * bvrs_other[j].s_blind0;// - alpha;


          norm_bvrs_b[j].blind     =   bvrs_b[j].b_blind ;// + bvr3[j].gamma[i];

          norm_bvrs_other[j].blind =   bvrs_other[j].b_blind;

          norm_bvrs_b[j].bb_u      =   bvrs_b[j].b_blind * bvrs_other[j].b_blind;// + bvr3[j].gamma[i];

          norm_bvrs_other[j].bb_u  =   bvrs_b[j].b_blind * bvrs_other[j].b_blind;// - bvr3[j].gamma[i];


          norm_progress[norm_step::blinds_created] = j + 1;
    }
}

    norm_bvrs * norm_bvrs_recv = (norm_bvrs *) std::aligned_alloc(sizeof(__m256i), nprofiles * sizeof(bvrs));

    void read_norm_bvrs(AES_KEY& key, tcp::socket& sin0, tcp::socket& sin1)
    {
      
      __m128i seed;

      read(sin0, boost::asio::buffer(&seed, sizeof(__m128i)));
     
      PRG(key, seed, (__m128i *)bvrs_created, (nprofiles * sizeof(bvrs)) / sizeof(__m128i));

      for(size_t j = 0; j < nprofiles; ++j)
      {
        read(sin1, boost::asio::buffer(&norm_bvrs_recv[j], sizeof(norm_bvrs_recv[j])));

        norm_progress[norm_step::norm_bvrs_in] = j + 1;
      }
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

        nitems_ = max_target;//dpfkey_read[j].nitems; 

        root_ = dpfkey_read[j].root;
         //std::cout << "-> root read: " << dpfkey_read[j].root.mX[0] << " " << dpfkey_read[j].root.mX[1] << std::endl; 
        for(size_t d = 0; d < depth; ++d) cw_[d] = dpfkey_read[j].cw[d];
        
        finalizer_ = dpfkey_read[j].finalizer;
        
        dpf_key<leaf_t, node_t, prgkey_t> dpfkey(std::move(nitems_),  std::move(root_), std::move(cw_), std::move(finalizer_), std::move(prgkey));

        keys_recv.push_back(dpfkey);

     
        norm_progress[norm_step::dpfs_in] = j + 1;
      }
    }

void write_blinded_sshare(size_t j, boost::asio::io_context &io_context, tcp::socket & sout)
{

	  while(norm_progress[norm_step::sshare_in] < j + 1 || norm_progress[norm_step::norm_bvrs_in] < j + 1) 
      { 
       std::this_thread::yield(); 
      }
     	
     blind_sshare[j]  = sshare_recv[j] + bvrs_created[j].s_blind1;// blind_recv[1];

      boost::asio::async_write(sout, boost::asio::buffer(&blind_sshare[j], sizeof(uint64_t)),
       [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      { 

        if(!ec)
        {
         if (j+1 < nprofiles) 
          {   
            write_blinded_sshare(j+1, io_context, sout); 
          }
        }
        else
        {
          write_blinded_sshare(j, io_context, sout);
        } 


    });
 
}

void read_blinded_sshare_recv(tcp::socket & sin)
{
  for(size_t j = 0; j < nprofiles; ++j)
  {	
    read(sin, boost::asio::buffer(&blinded_sshare_recv[j], sizeof(uint64_t)));
    norm_progress[norm_step::blind_sshare_in] = j + 1;
  }
}

void read_blinded_cnt_recv(tcp::socket & sin)
{
	for(size_t j = 0; j < nprofiles; ++j)
	{
  	 read(sin, boost::asio::buffer(&blinded_cnt_recv[j], sizeof(uint64_t)));
  	 norm_progress[norm_step::blind_cnt_in] = j + 1;
  	}
}

void write_z_to_P_b(size_t j, boost::asio::io_context &io_context, tcp::socket& sout)
{ 

	  while(norm_progress[norm_step::z_compute] < j + 1 ) 
    { 
       std::this_thread::yield(); 
    }

    io_context.restart();

    boost::asio::async_write(sout, boost::asio::buffer(&z[j], sizeof(int64_t)),
     [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
    { 

      if(!ec)
      {
       if (j+1 < nprofiles) 
        {   
          write_z_to_P_b(j+1, io_context, sout); 
        }
      }
      else
      {
        write_z_to_P_b(j, io_context, sout);
      } 


    });
  
  io_context.poll_one();
  norm_progress[norm_step::z_out] = j + 1;
	
}

void recv_shares_from_P_b(tcp::socket& sin)
{
   for(size_t j = 0; j < nprofiles; ++j)
   {
   	 read(sin, boost::asio::buffer(&X[j], sizeof(uint128_t)));
   	 norm_progress[norm_step::X_in] = j + 1;
   }
} 

void recv_reduced_from_P_b(tcp::socket& sin)
{
   for(size_t j = 0; j < nprofiles; ++j)
   {
     read(sin, boost::asio::buffer(&reduced[j], sizeof(uint64_t)));
     norm_progress[norm_step::reduced_in] = j + 1;
   }
} 

void send_s_to_P_b(size_t j, boost::asio::io_context& io_context, tcp::socket & sout)
{
 
   while(norm_progress[norm_step::s_compute] < j + 1)
   {
   	std::this_thread::yield();
   }
 
       
   async_write(sout, boost::asio::buffer(&s[j], sizeof(uint64_t)),
        [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
        
        {       
            if(!ec)
            {               
             if(j+ 1 < nprofiles) 
             {
             	send_s_to_P_b(j + 1, io_context, sout);
             }
      
            }          
            else
            {
              send_s_to_P_b(j , io_context, sout);
            }
        
        });
}

void write_blinded_cnt(size_t j, boost::asio::io_context &io_context, tcp::socket & sout)
{
 
  	while( norm_progress[norm_step::s_compute] < j + 1 || norm_progress[norm_step::norm_bvrs_in] < j + 1) 
     { 
       std::this_thread::yield(); 
     }

      blind_cnt[j]     = cnt[j] +  bvrs_created[j].s_blind0;// blind_recv[j];
    
   
      boost::asio::async_write(sout, boost::asio::buffer(&blind_cnt[j], sizeof(uint64_t)),
       [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      { 

        if(!ec)
        {
         if (j+1 < nprofiles) 
          {   
            write_blinded_cnt(j+1, io_context, sout); 
          }
        }
        else
        {
          write_blinded_cnt(j, io_context, sout);
        } 


    });

    norm_progress[norm_step::blind_cnt_out] = j + 1;
  
}

void read_sshare_recv(tcp::socket & sin)
{
  for(size_t j = 0; j < nprofiles; ++j)
  {	
   read(sin, boost::asio::buffer(&sshare_recv[j], sizeof(uint64_t)));
   norm_progress[norm_step::sshare_in] = j + 1;
  } 
}

void compute_z()
{
	for(size_t j = 0; j < nprofiles; ++j)
  	{	
  	   while(norm_progress[norm_step::s_compute] < j + 1 || norm_progress[norm_step::sshare_in] < j + 1 || 
  	   		 norm_progress[norm_step::blind_cnt_in] < j + 1 || norm_progress[norm_step::norm_bvrs_in] < j + 1 || norm_progress[norm_step::blind_sshare_in] < j + 1) 
      	{ 
      	  std::this_thread::yield(); 
      	}	
 	  
 	  z[j] = sshare_recv[j] * (cnt[j] + blinded_cnt_recv[j]) - bvrs_created[j].s_blind0 * (blinded_sshare_recv[j]) + norm_bvrs_recv[j].s_gamma;

 	  norm_progress[norm_step::z_compute] = j + 1;
	}
}