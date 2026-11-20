reset
set key right top
set xlabel "log T [K]"
set ylabel "log Fraction"
set xrange[1:6]
set yrange[-5:0]
set title "nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00"
p \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($17/$13)) w lp ti "e/H", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($18/$13)) w lp ti "HII/H", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($19/$13)) w lp ti "HI/H", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10(($13-$18-$19)/$13)) w lp ti "2*H2/H", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($20/$14)) w lp ti "HeIII/He", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($21/$14)) w lp ti "HeII/He", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($22/$14)) w lp ti "HeI/He", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($23/$15)) w lp ti "CIII/C", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($24/$15)) w lp ti "CII/C", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($25/$15)) w lp ti "CI/C", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($26/$16)) w lp ti "OIII/O", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($27/$16)) w lp ti "OII/O", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($28/$16)) w lp ti "OI/O", \
"out/PDRCoolingFunction_CO_nH1.0e+01_Nt1.0e+20_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10(($15-$23-$24-$25)/$15)) w lp ti "CO/C"
