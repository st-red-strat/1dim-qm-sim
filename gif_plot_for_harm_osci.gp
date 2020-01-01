set terminal gif animate optimize size 640, 480 enhanced
set output "./output/1dim-qm-sim.gif"
set xlabel "x"
set ylabel "|{/symbol y}|"
set samples 1000
filename = "./output/output.dat"

w=2*pi/300*(70.7*pi)

time_res = 20 # time resolution
tmax = 12800/10

H=0.004
dt=H**2/4
classic(int_t)=3.144660*(cos((60*pi*int_t*dt)-2.83709636))

do for[i=0:100] {
    set autoscale xy
    stats filename index 3*i*time_res using 1 nooutput
    set title sprintf("Potential Scattering, t = %04.0f", STATS_max )

    set yrange[-2:5]
    # set key bottom
 
    c=classic(STATS_max)
    set arrow 1 from c,-2 to c,5 nohead

    plot filename index 3*i*time_res+1 using 1:2  with lines title "Absolute value of wave function" ,\
    filename index 3*i*time_res+2 using 1:($2/w/1200) with lines title "Potential"

    unset arrow 1
}
