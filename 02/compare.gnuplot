#!/usr/bin/gnuplot --persist

data = ARG1
out = ARG2

set terminal png size 1920, 1080
set output out

plot data using 1:2 with lines title 'y_k' linewidth 4, \
    data using 1:3 with lines title 'y(x_k)' linewidth 4
