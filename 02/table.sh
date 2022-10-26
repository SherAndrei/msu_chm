#!/bin/bash

for num in 1 2 3 4 5 6; do
    for A in 1 10 1000; do
        echo -n "$num ";
        ./scheme$num.out $A;
    done;
done
