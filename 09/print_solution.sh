#!/bin/bash

die()
{
	echo >&2 error: "$@"
	exit 2
}

usage()
{
	cat >&2 <<EOF
Usage: ${0##*/} <input_data> <exact_function> <png>
DESCRIPTION:
  Construct png file with plotted exact solution and input data.
EOF
  exit 1
}

[ $# != 3 ] && usage

cat <<EOF | python3
from math import *;

with open("$1") as inp:
  with open("${1%%/*}/errors.txt", "w") as errors:
    for line in inp.readlines():
      x, y, val = map(float, line.split())
      error = abs(val - $2)
      errors.write(f"{x} {y} {error}\n")
EOF

cat <<EOF | gnuplot
set terminal png size 1600,800
set output '$3'

set hidden3d
set isosamples 30

set multiplot layout 1, 2

set label 1
  splot '$1' title 'Solution' lw 4, \
        $2 title 'Exact' lw 3
unset label 1

set label 2
  splot "${1%%/*}/errors.txt" title 'Errors' lw 2
unset label 2

unset multiplot
EOF
