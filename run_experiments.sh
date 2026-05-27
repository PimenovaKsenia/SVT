#!/bin/bash

echo "N,h,C_norm,L2_norm,Prec_time,Iter_time" > results.csv

for N in 20 40 80 120 160 200 248 512 1024
do
    echo "Запуск для N = $N ..."
    ./build/solver $N
done

python3 plot_results.py
