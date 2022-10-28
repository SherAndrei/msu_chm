#!/bin/bash

echo -n "№ "
echo -en "     E1     \t"
echo -en "     E2     \t"
echo -en "     E3     \t"
echo -en "     E6     \t"
echo "m A"
for scheme_num in 1 2 3 4 5 6; do
    for A in 1 10 1000; do
        echo -n "$scheme_num "
        for n in 1 2 3 6; do
            ./scheme$scheme_num.out $n $A;
            echo -en "\t"
        done;
        if [[ $scheme_num < 3 ]]; then
            echo -n "1 ";
        else
            echo -n "2 ";
        fi;
        echo $A
    done
done
