#ifndef P3_H__
#define P3_H__

#include "datatypes.h"
#include "network.h"
#include "aes.h"
#include "dpf.h"
#include <cmath> 
#include <random>


class server
{
public:
  server(boost::asio::io_context & io_context, const size_t nqueries, const AES_KEY * k)
    : nqueries_(nqueries),
      key(k),
      acceptor0(io_context, tcp::endpoint(tcp::v4(), PORT_P0_P3)),
      acceptor1(io_context, tcp::endpoint(tcp::v4(), PORT_P1_P3))
  {
    this->do_accept();
  }

private:
  friend class session;
  void do_accept();

  const size_t nqueries_;
  const AES_KEY * key;
  tcp::acceptor acceptor0;
  tcp::acceptor acceptor1;
};

class session : public std::enable_shared_from_this<session>
{
  public:

  session(const server * s);
  ~session();

  void start0(socket_ptr socket);
  void start1(socket_ptr socket);

  p3_beavers * p3_beavers_0;
  p3_beavers * p3_beavers_1;
  

  private:

  void write0(const size_t i, socket_ptr socket);
  void write1(const size_t i, socket_ptr socket);

  const size_t nqueries;

  const server * srvr;
 
};

#endif