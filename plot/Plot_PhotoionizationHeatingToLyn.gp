
reset
set logscale
set xrange[1:100]
set xlabel "Age [Myr]"
set ylabel "Heating Rate / N_Lyn"
p \
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($11/$20) w l ti "H I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($12/$21) w l ti "He I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($13/$22) w l ti "He II Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($14/$23) w l ti "C I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($15/$24) w l ti "C II Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($16/$25) w l ti "O I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($17/$26) w l ti "O II Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($18/$27) w l ti "Si I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($19/$28) w l ti "Si II Z = 0.0200"
