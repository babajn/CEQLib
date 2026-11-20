reset
set key left
set xlabel "log n [H/cc]"
set ylabel "log p/k [K/cc]"
set xrange[-3:+3]
set yrange[1:6]
p "out/PostProcessMap_Nt1.0e+19_ISRF0.3_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=0.3", "out/PostProcessMap_Nt1.0e+19_ISRF0.5_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=0.5", "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=1.0", "out/PostProcessMap_Nt1.0e+19_ISRF3.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=3.0", "out/PostProcessMap_Nt1.0e+19_ISRF10.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=10.0", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=100.0", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Nt=19 ISRF=1000.0"




