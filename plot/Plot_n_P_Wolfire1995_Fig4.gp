reset
set key left
set xlabel "log n [H/cc]"
set ylabel "log p/k [K/cc]"
set xrange[-3:+3]
set yrange[1:5]
p "out/PostProcessMap_Nt1.0e+18_ISRF1.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=18 ISRF=1.0", "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=1.0", "out/PostProcessMap_Nt1.0e+20_ISRF1.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=20 ISRF=1.0"
