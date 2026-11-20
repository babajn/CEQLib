reset
set key left bottom
set xlabel "log n [H/cc]"
set ylabel "log Fraction"
set xrange[-3:+3]
set yrange[-5:0]
p "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10($17/$13)) ti "e/H", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10($18/$13)) ti "Hii/H", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10($21/$14)) ti "Heii/He", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10($24/$15)) ti "Cii/C", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10($27/$16)) ti "Oii/O", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10(2*($13-$18-$19)/$13)) ti "2*H2/H", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00.dat" u (log10($13)):(log10(($15-$23-$24-$25)/$15)) ti "CO/C"
