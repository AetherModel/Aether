#!/bin/sh

g++ -I/usr/local/include -I/Users/ridley/Software/Json/json/include -o advect advect.cpp
./advect
./plot.py
open test.png
