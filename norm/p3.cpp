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



 AES_KEY prgkey;

 
 std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_read;
 std::vector<dpf_key<leaf_t, node_t, prgkey_t>> keys_recv;
 #include "isqrt.h"
 
 enum norm_step
 {
  blinds_created = 0,
  norm_bvrs_in,
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

         //target[j] = (X[j] >> ((P * precision)-6));
         target[j] =  mod(X[j] >> (3*precision-6), max_target); 
 

        // uint64_t * query_b = (uint64_t*)std::aligned_alloc(sizeof(__m128i), nitems / 8);
         query_b[j] = (uint64_t*)std::aligned_alloc(alignof(__m128i), keys_recv[j].full_bytes());
         
      
         keys_recv[j].evalfull((bool*)query_b[j]);
         
         auto [m, b, cnt0] = isqrt_coeffs(query_b[j], target[j]);
    
         cnt[j] = -cnt0;

         free(query_b[j]);
   
 
 
          
         long double int_, frac_ = std::modf(m, &int_);
 
         s[j] =  (uint64_t(-int_) * reduced[j]) + (uint64_t(-frac_*scale1_5) * reduced[j]) / scale1_5;
 
 
         norm_progress[norm_step::s_compute] = j + 1;
    }
}


 void keep_polling(boost::asio::io_context& io_context)
 {
    while(
             norm_progress[norm_step::norm_bvrs_in] < nprofiles || norm_progress[norm_step::norm_bvrs_P_other_out] < nprofiles 
          || norm_progress[norm_step::X_in] < nprofiles || norm_progress[norm_step::reduced_in] < nprofiles || norm_progress[norm_step::dpfs_in] < nprofiles 
          || norm_progress[norm_step::s_compute] < nprofiles || norm_progress[norm_step::sshare_in] < nprofiles || norm_progress[norm_step::blind_sshare_in] < nprofiles
          || norm_progress[norm_step::blind_cnt_in] < nprofiles ||  norm_progress[norm_step::z_compute] < nprofiles
          )
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
  
  boost::asio::io_context io_context;
    
  tcp::resolver resolver(io_context);
    
  const std::string host3 = (argc < 2) ? "127.0.0.1" : argv[1];



  tcp::acceptor acceptor1(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3));
  std::cout << " --- > " << std::endl;
  tcp::socket s1(acceptor1.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3 << std::endl;

  tcp::acceptor acceptor1_a(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_a));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_a(acceptor1_a.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_a << std::endl;

  tcp::acceptor acceptor1_b(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_b));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_b(acceptor1_b.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_b << std::endl;

  tcp::acceptor acceptor1_c(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_c));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_c(acceptor1_c.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_c << std::endl;

  tcp::acceptor acceptor1_d(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_d));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_d(acceptor1_d.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_d << std::endl;

  tcp::acceptor acceptor1_e(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_e));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_e(acceptor1_e.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_e << std::endl;

  tcp::acceptor acceptor1_f(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_f));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_f(acceptor1_f.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_f << std::endl;

  tcp::acceptor acceptor1_g(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_g));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_g(acceptor1_g.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_g << std::endl;

  tcp::acceptor acceptor1_i(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3_i));
  std::cout << " --- > " << std::endl;
  tcp::socket s1_i(acceptor1_i.accept());
  std::cerr << "Listenting on port: " << PORT_P1_P3_i << std::endl;

  tcp::acceptor acceptor0(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P3));
  std::cout << " --- > " << std::endl;
  tcp::socket s0(acceptor0.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P3 << std::endl;

  tcp::acceptor acceptor0_a(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P3_a));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_a(acceptor0_a.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P3_a << std::endl;

    tcp::acceptor acceptor0_b(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P3_b));
  std::cout << " --- > " << std::endl;
  tcp::socket s0_b(acceptor0_b.accept());
  std::cerr << "Listenting on port: " << PORT_P0_P3_b << std::endl;

  tcp::socket s2(io_context);
  boost::asio::connect(s2, resolver.resolve({host3,
  std::to_string(PORT_P3_P2)}));
  std::cout << "P0: [Established connection P3]" << std::endl;
  usleep(20000);

    tcp::socket s2_a(io_context);
  boost::asio::connect(s2_a, resolver.resolve({host3,
  std::to_string(PORT_P3_P2_a)}));
  std::cout << "P0: [Established connection P3]" << std::endl;
  usleep(20000);


    tcp::socket s2_b(io_context);
  boost::asio::connect(s2_b, resolver.resolve({host3,
  std::to_string(PORT_P3_P2_b)}));
  std::cout << "P0: [Established connection P3]" << std::endl;
  usleep(20000);

      tcp::socket s2_c(io_context);
  boost::asio::connect(s2_c, resolver.resolve({host3,
  std::to_string(PORT_P3_P2_c)}));
  std::cout << "P0: [Established connection P3]" << std::endl;
  usleep(20000);
 
      tcp::socket s2_d(io_context);
  boost::asio::connect(s2_d, resolver.resolve({host3,
  std::to_string(PORT_P3_P2_d)}));
  std::cout << "P0: [Established connection P3]" << std::endl;
  usleep(20000);
  try
  {

 
    std::thread poller(keep_polling, std::ref(io_context));

    std::thread blinds_creator(create_blinds, std::ref(key), std::ref(s0), std::ref(s2_a));
   
    std::thread blind_sender_to_Pb(send_blinds_to_P_b, 0, std::ref(io_context), std::ref(s0_b));

    std::thread blind_sender_to_P_other(send_blinds_to_P_other, 0, std::ref(io_context),  std::ref(s2_b));

    std::thread norm_reader(read_norm_bvrs, std::ref(key), std::ref(s2_d), std::ref(s2_c));

    std::thread share_from_Pb_receiver(recv_shares_from_P_b, std::ref(s1));

    std::thread reduced_from_Pb_receiver(recv_reduced_from_P_b, std::ref(s1_i));

    std::thread dpf_reader_from_Pb(read_dpfs_from_P_b, std::ref(s0_a));
   
    std::thread s_computer(compute_s);

    std::thread s_to_Pb_sender(send_s_to_P_b, 0, std::ref(io_context), std::ref(s1_a));
   
    std::thread sshare_receiver(read_sshare_recv, std::ref(s1_b));
 
    std::thread blinded_sshare_reader(read_blinded_sshare_recv, std::ref(s1_c));

    std::thread blinded_cnt_reader(read_blinded_cnt_recv, std::ref(s1_d));  

    std::thread blinded_sshare_writer(write_blinded_sshare, 0, std::ref(io_context), std::ref(s1_e)); 
   
    std::thread blinded_cnt_writer(write_blinded_cnt, 0, std::ref(io_context), std::ref(s1_f));  

    std::thread z_computer(compute_z);   
 
    std::thread z_to_P_b_writer(write_z_to_P_b, 0, std::ref(io_context), std::ref(s1_g));

    poller.join();
    
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

  
  }

  catch (std::exception & e)
  {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  return 0;
}
