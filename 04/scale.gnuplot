#!/usr/bin/gnuplot --persist

data = ARG1
out = ARG2

set terminal png size 1200,800
set output out

set title "Logarithmic dependence of the error on the number of steps"

set xlabel "log10(N)"
set ylabel "log10(1/||(u)_h-u_h||_h)"

set xrange [0:6]
set yrange [0:10]

plot data u (log10($1)):(log10(1/$2)) pointtype 5 notitle, \
    2*x with lines title 'y=2x'
