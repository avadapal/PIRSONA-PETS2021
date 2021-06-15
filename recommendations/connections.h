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

   #if (PARTY==0)
    // std::cerr << "P_ZERO: " << PORT_P0_P1 << std::endl;

    boost::asio::io_context io_context;
    
    tcp::resolver resolver(io_context);
    tcp::socket s2(io_context);
    boost::asio::connect(s2, resolver.resolve({host2,
    std::to_string(PORT_P0_P2)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);
    
    tcp::socket s2_a(io_context);
    boost::asio::connect(s2_a, resolver.resolve({host2,
    std::to_string(PORT_P0_P2_a)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);

    tcp::socket s2_b(io_context);
    boost::asio::connect(s2_b, resolver.resolve({host2,
    std::to_string(PORT_P0_P2_b)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);

    tcp::socket s2_c(io_context);
    boost::asio::connect(s2_c, resolver.resolve({host2,
    std::to_string(PORT_P0_P2_c)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);

    tcp::socket s2_d(io_context);
    boost::asio::connect(s2_d, resolver.resolve({host2,
    std::to_string(PORT_P0_P2_d)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);

    tcp::socket s2_e(io_context);
    boost::asio::connect(s2_e, resolver.resolve({host2,
    std::to_string(PORT_P0_P2_e)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);

    tcp::socket s2_f(io_context);
    boost::asio::connect(s2_f, resolver.resolve({host2,
    std::to_string(PORT_P0_P2_f)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;
    usleep(20000);

 
    tcp::socket sb(io_context);
    boost::asio::connect(sb, resolver.resolve({host1,
    std::to_string(PORT_P1_P0)}));
    // std::cerr << "P0: [Established connection P1 (1)]" << std::endl;
    usleep(20000);

    tcp::socket sb_a(io_context);
    boost::asio::connect(sb_a, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_a)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);
    
    tcp::socket sb_b(io_context);
    boost::asio::connect(sb_b, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_b)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);

    tcp::socket sb_c(io_context);
    boost::asio::connect(sb_c, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_c)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);

    tcp::socket sb_d(io_context);
    boost::asio::connect(sb_d, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_d)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);

    tcp::socket sb_e(io_context);
    boost::asio::connect(sb_e, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_e)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);

    tcp::socket sb_f(io_context);
    boost::asio::connect(sb_f, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_f)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);

    tcp::socket sb_g(io_context);
    boost::asio::connect(sb_g, resolver.resolve({host1,
    std::to_string(PORT_P1_P0_g)}));
    // std::cerr << "P0: [Established connection P1 (2)]" << std::endl;
    usleep(20000);
 
  
  #endif



  #if (PARTY==1)  
    // std::cerr << "P_ONE: " << PORT_P1_P0 << std::endl;
    boost::asio::io_context io_context;

    usleep(20000);
    tcp::resolver resolver(io_context);
    


    tcp::socket s2(io_context);
    boost::asio::connect(s2, resolver.resolve({host2,
    std::to_string(PORT_P1_P2)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
     usleep(20000); 

    tcp::socket s2_a(io_context);
    boost::asio::connect(s2_a, resolver.resolve({host2,
    std::to_string(PORT_P1_P2_a)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
    usleep(20000); 

    tcp::socket s2_b(io_context);
    boost::asio::connect(s2_b, resolver.resolve({host2,
    std::to_string(PORT_P1_P2_b)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
    usleep(20000); 
    
    tcp::socket s2_c(io_context);
    boost::asio::connect(s2_c, resolver.resolve({host2,
    std::to_string(PORT_P1_P2_c)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
    usleep(20000); 

    tcp::socket s2_d(io_context);
    boost::asio::connect(s2_d, resolver.resolve({host2,
    std::to_string(PORT_P1_P2_d)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
    usleep(20000); 

    tcp::socket s2_e(io_context);
    boost::asio::connect(s2_e, resolver.resolve({host2,
    std::to_string(PORT_P1_P2_e)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
    usleep(20000); 

    tcp::socket s2_f(io_context);
    boost::asio::connect(s2_f, resolver.resolve({host2,
    std::to_string(PORT_P1_P2_f)}));
    // std::cerr << "P1: [Established connection to P2]" << std::endl;    
    usleep(20000); 
 


    tcp::acceptor acceptor(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb(acceptor.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0 << std::endl;

    tcp::acceptor acceptor_a(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_a));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_a(acceptor_a.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_a << std::endl;
    

    tcp::acceptor acceptor_b(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_b));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_b(acceptor_b.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_b << std::endl;

    tcp::acceptor acceptor_c(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_c));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_c(acceptor_c.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_c << std::endl;

    tcp::acceptor acceptor_d(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_d));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_d(acceptor_d.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_d << std::endl;

    tcp::acceptor acceptor_e(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_e));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_e(acceptor_e.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_e << std::endl;

    tcp::acceptor acceptor_f(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_f));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_f(acceptor_f.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_f << std::endl;

    tcp::acceptor acceptor_g(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P0_g));
    // std::cerr << " --- > " << std::endl;
    tcp::socket sb_g(acceptor_g.accept());
    // std::cerr << "Listenting on port: " << PORT_P1_P0_g << std::endl;

   #endif
