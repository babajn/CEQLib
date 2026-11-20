reset
set key left
set xlabel "log n [H/cc]"
set ylabel "log p/k [K/cc]"
set xrange[-3:+3]
set yrange[1:5]
p \
    "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e-02.dat" u (log10($13)):(log10($9*$13)) ti "n-P Zs=0.01",\
    "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs3.0e-02.dat" u (log10($13)):(log10($9*$13)) ti "n-P Zs=0.03",\
    "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e-01.dat" u (log10($13)):(log10($9*$13)) ti "n-P Zs=0.1",\
    "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs3.0e-01.dat" u (log10($13)):(log10($9*$13)) ti "n-P Zs=0.3",\
    "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Zs=1",\
    "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs3.0e+00.dat" u (log10($13)):(log10($9*$13)) ti "n-P Zs=3"
