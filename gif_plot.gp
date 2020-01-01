set terminal gif animate optimize size 640, 480 enhanced
set output "./output/1dim-qm-sim.gif"
set xlabel "x"
set ylabel "|{/symbol y}|"
set samples 1000
filename = "./output/output.dat"

time_res = 5 # time resolution
tmax = 12800/10

do for[i=0:100] {
    set autoscale xy
    stats filename index 3*i*time_res using 1 nooutput
    set title sprintf("Potential Scattering, t = %04.0f", STATS_max )

    set yrange[-2:5]
    # set key bottom

    plot filename index 3*i*time_res+1 using 1:2  with lines title "Absolute value of wave function" ,\
    filename index 3*i*time_res+2 using 1:2 with lines title "Potential"
}
