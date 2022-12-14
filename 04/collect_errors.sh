#!/bin/bash

die()
{
	echo >&2 error: "$@"
	exit 2
}

usage() {
	cat >&2 <<EOF
Usage: ${0##*/} <Program Name>

Print dependendancy of ||(u)_h-u_h||_h from different N.
N is increasing exponentially by order of 10 starting from 10.
Errors is dumped to stdout in the next format
| N | ||(u)_h-u_h||_h |
EOF
	exit 2
}

check_prog()
{
	command -v "$1" >/dev/null 2>&1 || die "Program '$1' is required, but wasn't found."
}

[[ $# == 1 ]] || usage

check_prog awk
check_prog grep

PROG_PATH=$1

[[ -x "$PROG_PATH" ]] || die "$PROG_PATH does not exist or not an executable file"

for N in 10 100 1000 10000; do
    echo -en $N '\t';
    ./"$PROG_PATH" $N | grep 'Error' | awk '{print $2}';
done
