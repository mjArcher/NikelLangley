#!/bin/bash 
cd ./src/
g++ LibraryTest.cpp SquareTensor3.cpp ElasticPrimState.cpp -I ~/lib/ -I ~/lib/eigen/ -I ~/lib/eigen-dev/ -I ~/lib/eigen-dev/unsupported/ -lblitz -o test
./test

