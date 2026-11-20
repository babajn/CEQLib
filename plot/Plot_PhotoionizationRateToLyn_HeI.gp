
reset
set logscale
set xrange[1:100]
set xlabel "Age [Myr]"
set ylabel "Iozation Rate / N_Lyn"
set title "Photoionization Rate for He I"
p "out/Photoionization_spectra_Kroupa_Z0.0001.dat" u 1:($3/$21) w l ti "Z = 0.0001", "out/Photoionization_spectra_Kroupa_Z0.0004.dat" u 1:($3/$21) w l ti "Z = 0.0004", "out/Photoionization_spectra_Kroupa_Z0.0040.dat" u 1:($3/$21) w l ti "Z = 0.0040", "out/Photoionization_spectra_Kroupa_Z0.0080.dat" u 1:($3/$21) w l ti "Z = 0.0080", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:($3/$21) w l ti "Z = 0.0200", "out/Photoionization_spectra_Kroupa_Z0.0500.dat" u 1:($3/$21) w l ti "Z = 0.0500", "out/Photoionization_spectra_Kroupa_Z0.1000.dat" u 1:($3/$21) w l ti "Z = 0.1000"
