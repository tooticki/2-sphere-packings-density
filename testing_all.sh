#!/bin/bash

./test 1 1 1 1 1 > 1111.txt

./test 1 1 r r 1 > 11rr.txt
./test 1 r 1 r 1 > 1r1r.txt
./test r r 1 1 1 > rr11.txt

./test 1 1 1 r 1 > 111r.txt
./test 1 r 1 1 1 > 1r11.txt
