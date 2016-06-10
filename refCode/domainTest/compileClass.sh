#! /bin/bash
g++ -c ~/solid-2D/src/ElasticPrimState.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz &> ~/solid-2D/domainTest/compile.out
g++ -c ~/solid-2D/src/ElasticState.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz &> ~/solid-2D/domainTest/compile.out
g++ -c ~/solid-2D/src/SquareTensor3.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz &> ~/solid-2D/domainTest/compile.out
g++ -c ~/solid-2D/src/Domain.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz &> ~/solid-2D/domainTest/compile.out
g++ -o domainTest Domain.o SquareTensor3.o ElasticState.o ElasticPrimState.o -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz &> ~/solid-2D/domainTest/compile.out

./domainTest
