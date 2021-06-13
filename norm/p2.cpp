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
 
#include<fstream>
#include "aes.h"
#include "network.h"
 

 #include "dpf2.h"

 using namespace dpf;

typedef bool leaf_t;
typedef block<__m128i> node_t;
typedef AES_KEY prgkey_t;



using boost::asio::ip::tcp;

boost::asio::io_context io_context;
 
 prgkey_t prgkey;

   std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_read;
   std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_recv;
 #include "isqrt.h"
 
 
 
  uint64_t key0 = 597349;
  uint64_t key1 = 121379;

 enum norm_step
 {
  blinds_created = 0,
  norm_bvrs_in,
  dotprod_gen,
  dotprod_bvr0_out,
  dotprod_bvr1_out,
  X_in,
  reduced_in,
  norm_bvrs_PB_out,
  norm_bvrs_P_other_out,
  dpfs_in,
  sshare_in,
  blind_cnt_in, 
  blind_sshare_in,
  blind_cnt_out,
  s_compute,
  z_compute,
  z_out,
  norm_num_steps,
};
size_t norm_progress[norm_step::norm_num_steps] = { 0 };
#include "p_other.h"


uint64_t * query_b[nprofiles];


 void compute_s()
 {
  
    for(size_t j = 0; j < nprofiles; ++j)
    {
      
         while(norm_progress[norm_step::dpfs_in] < j + 1 || norm_progress[norm_step::X_in] < j + 1 || norm_progress[norm_step::reduced_in] < j + 1) 
         { 
 
            std::this_thread::yield(); 
         }
 
         target[j] =  mod(X[j] >> (3*precision-6), max_target); 
  
         query_b[j] = (uint64_t*)std::aligned_alloc(alignof(__m128i), keys_recv[j].full_bytes());
         keys_recv[j].evalfull((bool*)query_b[j]);

         auto [m, b, cnt0] = isqrt_coeffs(query_b[j], target[j]);
 
         cnt[j] = -cnt0;

         free(query_b[j]);
 
         long double int_, frac_ = std::modf(m, &int_);
 
         s[j] = uint64_t(b) + (uint64_t(-int_) * reduced[j]) + (uint64_t(-frac_*scale1_5) * reduced[j]) / scale1_5;
 
         norm_progress[norm_step::s_compute] = j + 1;
    } 
 }
void create_norm_beavers_usr()
{
  for(size_t j = 0; j < nprofiles; ++j)
  {
     arc4random_buf(&bvr3[j], sizeof(bvr3[j]));

     dot_bvr0[j].u_blind = bvr3[j].u0_blind;
     dot_bvr0[j].s       = dot(bvr3[j].u0_blind, bvr3[j].u1_blind) + bvr3[j].alpha;;
 
     dot_bvr1[j].u_blind = bvr3[j].u1_blind;
     dot_bvr1[j].s       = dot(bvr3[j].u0_blind, bvr3[j].u1_blind) - bvr3[j].alpha;
     
     dot_bvr0[j]._128_u_blind = bvr3[j]._128_u0_blind;
     dot_bvr0[j]._128_s       = dot(bvr3[j]._128_u0_blind, bvr3[j]._128_u1_blind) + bvr3[j]._128_alpha;;
 
     dot_bvr1[j]._128_u_blind = bvr3[j]._128_u1_blind;
     dot_bvr1[j]._128_s       = dot(bvr3[j]._128_u0_blind, bvr3[j]._128_u1_blind) - bvr3[j]._128_alpha;

     norm_progress[norm_step::dotprod_gen] = j + 1;
  }

  
}

void write_dotprod_bvrs_to_P0(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{
      while(norm_progress[norm_step::dotprod_gen] < j + 1)
    {
      std::this_thread::yield();
    } 

      boost::asio::async_write(sout, boost::asio::buffer(&dot_bvr0[j], dotprod_stuff::size),
      [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      { 

        if(!ec)
        {
         if (j+1 < nprofiles) 
          {   
            write_dotprod_bvrs_to_P0(j+1, io_context, sout); 
          }
        }
        else
        {
          write_dotprod_bvrs_to_P0(j, io_context, sout);
        } 


    });

      norm_progress[norm_step::dotprod_bvr0_out] = j + 1;

}


void write_dotprod_bvrs_to_P1(size_t j, boost::asio::io_context& io_context, tcp::socket& sout)
{

    while(norm_progress[norm_step::dotprod_gen] < j + 1)
    {
      std::this_thread::yield();
    }
 
      boost::asio::async_write(sout, boost::asio::buffer(&dot_bvr1[j], dotprod_stuff::size),
       [j, &io_context, &sout](boost::system::error_code ec, std::size_t) 
      { 

        if(!ec)
        {
         if (j+1 < nprofiles) 
          {   
            write_dotprod_bvrs_to_P1(j+1, io_context, sout); 
          }
        }
        else
        {
          write_dotprod_bvrs_to_P1(j, io_context, sout);
        } 


    });

      norm_progress[norm_step::dotprod_bvr1_out] = j + 1;


}

 void keep_polling(boost::asio::io_context& io_context)
 {
    while(norm_progress[norm_step::norm_bvrs_in]     < nprofiles || norm_progress[norm_step::dotprod_bvr0_out] < nprofiles || 
          norm_progress[norm_step::dotprod_bvr1_out] < nprofiles || norm_progress[norm_step::norm_bvrs_P_other_out] < nprofiles || norm_progress[norm_step::X_in] < nprofiles ||
          norm_progress[norm_step::reduced_in] < nprofiles || norm_progress[norm_step::dpfs_in] < nprofiles 
          || norm_progress[norm_step::s_compute] < nprofiles || norm_progress[norm_step::sshare_in] < nprofiles || norm_progress[norm_step::blind_sshare_in] < nprofiles
          || norm_progress[norm_step::blind_cnt_in] < nprofiles || norm_progress[norm_step::z_compute] < nprofiles)
    {
      io_context.reset();
      io_context.poll();
    }
 }


int main(int argc, char* argv[])
{ 
  
  keys_read.reserve(nprofiles);
  keys_recv.reserve(nprofiles);
  
  AES_KEY key;
 

 

  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &key);
  
  tcp::acceptor acceptor1(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2));
  std::cout << " --- > " << std::endl;
  tcp::socket s1(acceptor1.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P2 << std::endl;

  tcp::acceptor acceptor1_a(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_a));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_a(acceptor1_a.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P2_a << std::endl;

  tcp::acceptor acceptor1_b(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_b));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_b(acceptor1_b.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P2_b << std::endl;


  tcp::acceptor acceptor1_h(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P2_h));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_dotprod(acceptor1_h.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P2_h << std::endl;

  tcp::acceptor acceptor0(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2));
  std::cout << " --- > " << std::endl;
  tcp::socket s0(acceptor0.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2 << std::endl;

  tcp::acceptor acceptor0_a(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_a));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_a(acceptor0_a.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_a << std::endl;

  tcp::acceptor acceptor0_b(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_b));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_b(acceptor0_b.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_b << std::endl;

  tcp::acceptor acceptor0_c(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_c));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_c(acceptor0_c.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_c << std::endl;

  tcp::acceptor acceptor0_d(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_d));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_d(acceptor0_d.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_d << std::endl;

  tcp::acceptor acceptor0_e(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_e));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_e(acceptor0_e.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_e << std::endl;

  tcp::acceptor acceptor0_f(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_f));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_f(acceptor0_f.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_f << std::endl;

  tcp::acceptor acceptor0_g(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_g));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_g(acceptor0_g.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_g << std::endl;

    tcp::acceptor acceptor0_i(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_i));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_i(acceptor0_i.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_i << std::endl;

  tcp::acceptor acceptor0_h(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P2_h));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_dotprod(acceptor0_h.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P2_h << std::endl;
  
  tcp::acceptor acceptor3(io_context, tcp::endpoint(tcp::v4(), PORT_P3_P2));
  std::cout << " --- > " << std::endl;
  tcp::socket s3(acceptor3.accept());
  std::cerr << "Listenting on port: " << PORT_P3_P2 << std::endl;

  tcp::acceptor acceptor3_a(io_context, tcp::endpoint(tcp::v4(), PORT_P3_P2_a));
  std::cout << " --- > " << std::endl;
  tcp::socket s3_a(acceptor3_a.accept());
  std::cerr << "Listenting on port: " << PORT_P3_P2_a << std::endl;

    tcp::acceptor acceptor3_b(io_context, tcp::endpoint(tcp::v4(), PORT_P3_P2_b));
  std::cout << " --- > " << std::endl;
  tcp::socket s3_b(acceptor3_b.accept());
  std::cerr << "Listenting on port: " << PORT_P3_P2_b << std::endl;
 
     tcp::acceptor acceptor3_c(io_context, tcp::endpoint(tcp::v4(), PORT_P3_P2_c));
  std::cout << " --- > " << std::endl;
  tcp::socket s3_c(acceptor3_c.accept());
  std::cerr << "Listenting on port: " << PORT_P3_P2_c << std::endl;


  tcp::acceptor acceptor3_d(io_context, tcp::endpoint(tcp::v4(), PORT_P3_P2_d));
  std::cout << " --- > " << std::endl;
  tcp::socket s3_d(acceptor3_d.accept());
  std::cerr << "Listenting on port: " << PORT_P3_P2_d << std::endl;

   std::thread poller(keep_polling, std::ref(io_context));

   std::thread dotprod_bvr_creater(create_norm_beavers_usr);

   std::thread dotprod_to_P0_writer(write_dotprod_bvrs_to_P0, 0, std::ref(io_context), std::ref(s0_dotprod));

   std::thread dotprod_to_P1_writer(write_dotprod_bvrs_to_P1, 0, std::ref(io_context), std::ref(s1_dotprod));

   std::thread blinds_creator(create_blinds, std::ref(key), std::ref(s1), std::ref(s3_d));
   
   std::thread blind_sender_to_Pb(send_blinds_to_P_b, 0, std::ref(io_context),  std::ref(s1_b));

   std::thread blind_sender_to_P_other(send_blinds_to_P_other, 0, std::ref(io_context), std::ref(s3_c));

   std::thread norm_reader(read_norm_bvrs, std::ref(key), std::ref(s3_a), std::ref(s3_b));
   
   std::thread share_from_Pb_receiver(recv_shares_from_P_b, std::ref(s0));

   std::thread reduced_from_Pb_receiver(recv_reduced_from_P_b, std::ref(s0_i));

    std::thread dpf_reader_from_Pb(read_dpfs_from_P_b, std::ref(s1_a));

   std::thread s_computer(compute_s);
   
    std::thread s_to_Pb_sender(send_s_to_P_b, 0, std::ref(io_context), std::ref(s0_a));
   
    std::thread sshare_receiver(read_sshare_recv, std::ref(s0_b));

    std::thread blinded_sshare_reader(read_blinded_sshare_recv, std::ref(s0_c));

    std::thread blinded_cnt_reader(read_blinded_cnt_recv, std::ref(s0_d));  

   std::thread blinded_sshare_writer(write_blinded_sshare, 0, std::ref(io_context), std::ref(s0_e)); 
   
   std::thread blinded_cnt_writer(write_blinded_cnt, 0, std::ref(io_context), std::ref(s0_f));  

   std::thread z_computer(compute_z);   
   
   std::thread z_to_P_b_writer(write_z_to_P_b,  0, std::ref(io_context), std::ref(s0_g));
  

   poller.join();

   dotprod_bvr_creater.join();
   
   dotprod_to_P0_writer.join();
  
   dotprod_to_P1_writer.join();

   blinds_creator.join();
   
   blind_sender_to_Pb.join();
   
   blind_sender_to_P_other.join();  
   
   norm_reader.join();
   
   share_from_Pb_receiver.join();
   
   reduced_from_Pb_receiver.join();

   dpf_reader_from_Pb.join();
  
   s_computer.join();
   
   s_to_Pb_sender.join();

   sshare_receiver.join();

   blinded_sshare_reader.join();
   
   blinded_cnt_reader.join();

   blinded_sshare_writer.join();
   
   blinded_cnt_writer.join();
   
   z_computer.join();

   z_to_P_b_writer.join();

 
  
  return 0;
}

 