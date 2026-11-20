reset
set logscale
#set xrange[1:100]
set xlabel "Z"
set ylabel "Ionization Rate / N_Lyn"
p \
"out/PhotoionizationNormalizedNLyn.dat" u 1:3 w l ti "H I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:6 w l ti "He I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:9 w l ti "He II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:12 w l ti "C I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:15 w l ti "C II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:18 w l ti "O I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:21 w l ti "O II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:24 w l ti "Si I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:27 w l ti "Si II"
