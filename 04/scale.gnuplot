#!/usr/bin/gnuplot --persist

data = ARG1

set xlabel "log10(N)"
set ylabel "log10(1/||(u)_h-u_h||_h)"

set xrange [0:6]
set yrange [0:10]

plot data u (log10($1)):(log10(1/$2)) pointtype 5 notitle, \
    2*x with lines title 'y=2x'

