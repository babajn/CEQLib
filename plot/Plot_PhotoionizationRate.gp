
reset
set logscale
PC_CGS = 3.08568025e+18
d = 10.0*PC_CGS
fact = 1.0/4.0/3.14/d**2
set xrange[1:300]
set xlabel "Age [Myr]"
set ylabel "[1/s/Msun] for d = 10 pc"
p "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$2) w l ti "H I", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$3) w l ti "He I", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$4) w l ti "He II", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$5) w l ti "C I", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$6) w l ti "C II", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$7) w l ti "O I", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$8) w l ti "O II", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$9) w l ti "Si I", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$10) w l ti "Si II"
