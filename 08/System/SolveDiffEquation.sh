#!/usr/bin/bash

die()
{
	echo >&2 error: "$@"
	exit 2
}

usage() {
	cat >&2 <<EOF
Usage: ${0##*/} [-e <1e-10>] [-m <10>] [-X <1.0>] [-a <0>] [-b <0>]

DESCRIPTION:
  Plot solution of differential equation y''=f(x,y)
  for x in (0, X), y(0) = a, y(X) = b, found by
  solving system of non-linear algebraic equations
  with Newton's method.

OPTIONS:
  -h                -- see help
  -e ( = 1e-10    ) -- precision of the solution
  -m ( = 10       ) -- amounf of equations in the system, m > 0
  -X ( = 1.0      ) -- right bound for x
  -a ( = 0        ) -- value of y(0)
  -b ( = 0        ) -- value of y(X)
  -o ( = 'out.png') -- png_output file with solution
EOF
	exit 2
}

check_prog()
{
	command -v "$1" >/dev/null 2>&1 || die "Program '$1' is required, but wasn't found."
}

check_file()
{
  [ -f $1 ] || die "File '$1' is expected, but wasn't found."
}

check_arg_for_regex()
{
  ([ ! -z $1 ] && [[ $1 =~ $2 ]]) || die "input argument didn't match the requirements: '$1'"
}

convert_to_txt_output()
{
  local output_dir=$(dirname $1)
  local dirless_name=${1##*/}
  local pure_output_name=${dirless_name%.*}
  local txt_output="$output_dir/${pure_output_name}_solution.txt"
  echo "$txt_output"
}

strict_non_zero_unsigned='^\+?[1-9][0-9]*$'

double_impl='[0-9]+(\.(([0-9]+)?))?((e|E)((\+|-)?)[0-9]+)?'
strict_double="^(\+|-)?$double_impl$"
strict_unsigned_double="^\+?$double_impl$"

dir=$(dirname $0)
check_file "${dir}/DiffEquation.c"
check_file "${dir}/Newton.c"
check_file "${dir}/Makefile"

check_prog paste
check_prog bc
check_prog sed
check_prog make
check_prog gnuplot

eps=1e-10
m=10
X=1.0
a=0.0
b=0.0
png_output='out.png'

while getopts "a:b:e:m:X:o:h" flag
do
    case "${flag}" in
        e) eps=${OPTARG}
          check_arg_for_regex $eps $strict_unsigned_double
          ;;
        m) m=${OPTARG}
          check_arg_for_regex $m $strict_non_zero_unsigned
          ;;
        X) X=${OPTARG}
          check_arg_for_regex $X $strict_unsigned_double
          ;;
        a) a=${OPTARG}
          check_arg_for_regex $a $strict_double
          ;;
        b) b=${OPTARG}
          check_arg_for_regex $a $strict_double
          ;;
        o) png_output=${OPTARG};;
        h|*) usage;;
    esac
done

sed -i -E "s/double X = .+?;/double X = $X;/g" "${dir}/DiffEquation.c" || die "cannot match X"
sed -i -E "s/double a = .+?;/double a = $a;/g" "${dir}/DiffEquation.c" || die "cannot match a"
sed -i -E "s/double b = .+?;/double b = $b;/g" "${dir}/DiffEquation.c" || die "cannot match b"
make -C "$dir" > /dev/null

newton_res=$("$dir"/DiffEquation.out $eps $m)
newton_error_code=$?

if [ $newton_error_code != 0 ] ; then
  die "Newton failed with error: " $newton_res ", code (" $newton_error_code ")";
fi

h=$(bc -l <<< "scale=14;$X/$(($m + 1))")
last=$(bc -l <<< "scale=14;$X - $h")

# with initial conditions
array_of_x="0\\n$(seq $h $h $last)\\n$X"
array_of_y="$a\\n$(echo -e "$newton_res")\\n$b"

txt_output=$(convert_to_txt_output $png_output)
paste <(echo -e "$array_of_x") <(echo -e "$array_of_y") > "$txt_output"

gnuplot <<EOF
set terminal png size 1200,800
set output "$png_output"

plot "$txt_output" with lp title 'Solution' lw 4
EOF
