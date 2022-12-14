#!/usr/bin/gnuplot --persist

data = ARG1
out = ARG2

set terminal png size 1000,800
set output out

set xrange [0:1]

plot data using 1:2 with linespoints title 'Method', \
    data using 1:3 with linespoints title 'Original'
