#! /bin/bash
g++ -c ~/solid-2D/domainTest/src/ElasticPrimState.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz 
g++ -c ~/solid-2D/domainTest/src/SquareTensor3.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz 
g++ -c ~/solid-2D/domainTest/src/DomainEigen.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz 
g++ -o domainTest DomainEigen.o SquareTensor3.o ElasticPrimState.o -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz 
# g++ ./src/DomainEigen.cpp ~/solid-2D/src/SquareTensor3.cpp ~/solid-2D/src/ElasticPrimState.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz -o domainTest
./domainTest
