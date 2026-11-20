reset
set logscale
#set xrange[1:100]
set xlabel "Z"
set ylabel "Photons / N_Lyn"
p \
"out/PhotoionizationNormalizedNLyn.dat" u 1:2 w l ti "H I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:5 w l ti "He I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:8 w l ti "He II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:11 w l ti "C I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:14 w l ti "C II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:17 w l ti "O I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:20 w l ti "O II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:23 w l ti "Si I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:26 w l ti "Si II"
