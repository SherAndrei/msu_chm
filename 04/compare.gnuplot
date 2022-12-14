#!/usr/bin/gnuplot --persist

data = ARG1

set xrange [0:1]

plot data using 1:2 title 'Method', \
    data using 1:3 title 'Original'
