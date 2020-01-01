#!/bin/sh

g++ -std=c++14 -Wall -o 1dim-qm-sim ./1dim-qm-sim.cpp -lm
./1dim-qm-sim 1>output/output.dat 2>output/error.dat
# cat output/output.dat
# cat output/error.dat
gnuplot gif_plot.gp
# gnuplot gif_plot_for_harm_osci.gp
