
reset
set logscale
PC_CGS = 3.08568025e+18
d = 10.0*PC_CGS
fact = 1.0/4.0/3.14/d**2
set xrange[1:1000]
set xlabel "Age [Myr]"
set ylabel "Heating Rate [erg*cm3/s/Msun] for d = 10 pc"
set title "Photoionization Heating Rate for Hi"
p "out/Photoionization_spectra_Kroupa_Z0.0001.dat" u 1:(fact*$11) w l ti "Z = 0.0001", "out/Photoionization_spectra_Kroupa_Z0.0004.dat" u 1:(fact*$11) w l ti "Z = 0.0004", "out/Photoionization_spectra_Kroupa_Z0.0040.dat" u 1:(fact*$11) w l ti "Z = 0.0040", "out/Photoionization_spectra_Kroupa_Z0.0080.dat" u 1:(fact*$11) w l ti "Z = 0.0080", "out/Photoionization_spectra_Kroupa_Z0.0200.dat" u 1:(fact*$11) w l ti "Z = 0.0200", "out/Photoionization_spectra_Kroupa_Z0.0500.dat" u 1:(fact*$11) w l ti "Z = 0.0500", "out/Photoionization_spectra_Kroupa_Z0.1000.dat" u 1:(fact*$11) w l ti "Z = 0.1000"
