#! /bin/bash
g++ Input.cpp -I ~/lib/eigen/  -I ~/lib/libconfig/lib/ -lconfig++ -o configTest
# -L /home/raid/ma595/install/Libraries/config4cpp/lib/ -lconfig4cpp -I /home/raid/ma595/install/Libraries/config4cpp/include/ 
g++ test.cpp -I ~/lib/eigen/  -I ~/lib/libconfig/lib/ -lconfig++ -o configTest
./configTest
