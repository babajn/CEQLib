reset
set key left
set xlabel "log n [H/cc]"
set ylabel "log Cooling/Heating Rates [erg/s/H]"
set xrange[-3:+3]
set yrange[-28:-22]
p "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($33/$13)) ti "L_[Cii]158", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($35/$13)) ti "L_[Oi ]063", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($36/$13)) ti "L_DustRec", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($37/$13)) ti "L_H2", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($38/$13)) ti "L_CO", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($39/$13)) ti "H_PE", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($40/$13)) ti "H_CR", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($41/$13)) ti "H_XR", "out/PostProcessMap_Nt1.0e+19_ISRF1000.0_Zs1.0e+00.dat" u (log10($13)):(log10($42/$13)) ti "H_H2"
