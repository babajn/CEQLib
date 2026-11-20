reset
set key left
set xlabel "log n/X0 [H/cc]"
set ylabel "log P/k/X0 [K/cc]"
set xrange[-3:+3]
set yrange[1:5]
p "out/PostProcessMap_Nt1.0e+19_ISRF0.1_Zs1.0e+00_Xo0.1.dat" u (log10($13/0.1)):(log10($9*$13/0.1)) ti "n-P X0=0.1", "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e+00_Xo1.0.dat" u (log10($13/1.0)):(log10($9*$13/1.0)) ti "n-P X0=1.0", "out/PostProcessMap_Nt1.0e+19_ISRF10.0_Zs1.0e+00_Xo10.0.dat" u (log10($13/10.0)):(log10($9*$13/10.0)) ti "n-P X0=10.0", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00_Xo100.0.dat" u (log10($13/100.0)):(log10($9*$13/100.0)) ti "n-P X0=100.0", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00_Xo1000.0.dat" u (log10($13/1000.0)):(log10($9*$13/1000.0)) ti "n-P X0=1000.0"
#p "out/PostProcessMap_Nt1.0e+19_ISRF0.1_Zs1.0e+00_Xo0.1.dat" u (log10($13)):(log10($9*$13)) ti "n-P X0=0.1", "out/PostProcessMap_Nt1.0e+19_ISRF1.0_Zs1.0e+00_Xo1.0.dat" u (log10($13)):(log10($9*$13)) ti "n-P X0=1.0", "out/PostProcessMap_Nt1.0e+19_ISRF10.0_Zs1.0e+00_Xo10.0.dat" u (log10($13)):(log10($9*$13)) ti "n-P X0=10.0", "out/PostProcessMap_Nt1.0e+19_ISRF100.0_Zs1.0e+00_Xo100.0.dat" u (log10($13)):(log10($9*$13)) ti "n-P X0=100.0", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00_Xo1000.0.dat" u (log10($13)):(log10($9*$13)) ti "n-P X0=1000.0"
