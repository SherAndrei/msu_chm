#!/usr/bin/gnuplot --persist

data = ARG1
out = ARG2

set terminal png size 1200,800
set output out

plot data u 1:2 w lp, data u 1:3 w lp
