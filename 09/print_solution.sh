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

error_filename="${1%%.*}_errors.txt"

min_max_mean=$(cat <<EOF | python3
from math import *;

with open("$1") as inp:
  with open("$error_filename", "w") as errors:
    all_errors = []
    for line in inp.readlines():
      x, y, val = map(float, line.split())
      error = abs($2 - val)
      all_errors.append(error)
      errors.write(f"{x} {y} {error}\n")
    print("min:", min(all_errors), "max:", max(all_errors), "mean:", sum(all_errors) / len(all_errors))
EOF
)

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
  set title "$min_max_mean"
  splot "$error_filename" title 'Errors' lw 2
unset label 2

unset multiplot
EOF
