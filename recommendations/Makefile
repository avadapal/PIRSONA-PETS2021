#
# Copyright (C) 2019  Anonymous
#
# This is a pre-release version of the DPF++ library distributed anonymously
# for peer review. A public release of the software will be published under the
# LPGL v2.1 license in the near future. Please do not redistribute this version
# of the software.
#
CXX=g++
CXXFLAGS=--std=c++17 -march=native -O3 -g -Wno-ignored-attributes -fopenmp -DBOOST_ERROR_CODE_HEADER_ONLY #-DDEBUG #-DVERBOSE
LIBS=-lbsd -lcrypto -lssl -lboost_system -lboost_thread -pthread

default: p0 p1 p2 

all: parties

p0: pB.cpp dpf.h prg.h compare.h connections.h zero-out.h dotproduct.h
	$(CXX) $(CXXFLAGS) -o p0 pB.cpp $(LIBS) -DPARTY=0 -I.

p1: pB.cpp dpf.h prg.h compare.h connections.h zero-out.h dotproduct.h
	$(CXX) $(CXXFLAGS) -o p1 pB.cpp $(LIBS) -DPARTY=1 -I.

p2: p2.cpp network.h prg.h aes.h compare.h connections.h zero-out.h dotproduct.h
	$(CXX) $(CXXFLAGS) -o p2 p2.cpp $(LIBS) -I.

server: server.cpp network.h prg.h aes.h connections.h
	$(CXX) $(CXXFLAGS) -o server server.cpp $(LIBS) -I.

client: client.cpp network.h prg.h aes.h connections.h
	$(CXX) $(CXXFLAGS) -o client client.cpp $(LIBS) -I.

parities: parities.cpp dpf.h prg.h block.h
	g++ -ggdb -O3 -march=native -std=c++17 -Wall -pedantic -o parities parities.cpp -lbsd -Wno-ignored-attributes

lowmc.o: lowmc.cpp lowmc.h block.h
	g++ -ggdb -O3 -march=native -std=c++17 -Wall -pedantic -c lowmc.cpp -lbsd -Wno-ignored-attributes

clean:
	rm -f p0 p1 p2 p3

kill:
	pkill p2 & pkill p1 & pkill p0 & pkill p3
