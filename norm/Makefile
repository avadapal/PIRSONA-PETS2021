CXX=g++
CXXFLAGS=--std=c++17 -march=native -O3 -g -Wno-ignored-attributes -fopenmp -DBOOST_ERROR_CODE_HEADER_ONLY #-DDEBUG #-DVERBOSE
LIBS=-lbsd -lboost_system -lboost_thread -pthread

default: p0 p1 p2 p3 

p0: pB.cpp dpf.h prg.h isqrt.h connections.h
	$(CXX) $(CXXFLAGS) -o p0 pB.cpp $(LIBS) -DPARTY=0 -I.

p1: pB.cpp dpf.h prg.h isqrt.h connections.h
	$(CXX) $(CXXFLAGS) -o p1 pB.cpp $(LIBS) -DPARTY=1 -I.

p2: p2.cpp network.h prg.h aes.h fixed.h isqrt.h p_other.h connections.h
	$(CXX) $(CXXFLAGS) -o p2 p2.cpp $(LIBS) -I.

p3: p3.cpp network.h prg.h aes.h fixed.h isqrt.h p_other.h connections.h
	$(CXX) $(CXXFLAGS) -o p3 p3.cpp $(LIBS) -I.

isqrt2: isqrt2.cpp dpf.h prg.h 
	$(CXX) $(CXXFLAGS) -o isqrt2 isqrt2.cpp $(LIBS) -DPARTY=0 -I.

gen_users: gen_users.cpp datatypes.h fixed.h common.h
	touch uprofiles
	touch uprofilesA
	touch uprofiles.0
	touch uprofiles.1
	$(CXX) $(CXXFLAGS) -o gen_users gen_users.cpp $(LIBS)

clean:
	rm -f p0 p1 p2 p3  

kill:
	pkill p2 & pkill p1 & pkill p0 & pkill p3

