
reset
set logscale
set xrange[1:100]
set xlabel "Age [Myr]"
set ylabel "Iozation Rate / N_Lyn"
p \
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($2/$20) w l ti "H I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($3/$21) w l ti "He IZ = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($4/$22) w l ti "He II Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($5/$23) w l ti "C I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($6/$24) w l ti "C II Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($7/$25) w l ti "O I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($8/$26) w l ti "O II Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($9/$27) w l ti "Si I Z = 0.0200",\
"out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($10/$28) w l ti "Si II Z = 0.0200"



