reset
set logscale
#set xrange[1:100]
set xlabel "Z"
set ylabel "Heating Rate / N_Lyn"
p \
"out/PhotoionizationNormalizedNLyn.dat" u 1:4 w l ti "H I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:7 w l ti "He I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:10 w l ti "He II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:13 w l ti "C I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:16 w l ti "C II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:19 w l ti "O I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:22 w l ti "O II",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:25 w l ti "Si I",\
"out/PhotoionizationNormalizedNLyn.dat" u 1:28 w l ti "Si II"
