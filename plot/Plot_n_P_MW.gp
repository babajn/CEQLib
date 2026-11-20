reset
set key left
set xlabel "log n [H/cc]"
set ylabel "log p/k [K/cc]"
set xrange[-3:+3]
set yrange[0:6]
p "out/PostProcessMap_Nt1.0e+18_Go1.7e+00_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=18 G0=1.7", "out/PostProcessMap_Nt1.0e+19_Go1.7e+00_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 G0=1.7", "out/PostProcessMap_Nt1.0e+20_Go1.7e+00_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=20 G0=1.7", "out/PostProcessMap_Nt1.0e+20_Go1.7e+01_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=20 G0=17."
