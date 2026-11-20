reset
set key left
set xlabel "log T [K]"
set ylabel "log Cooling,Heating/nH^2 [erg cm3/s]"
set xrange[1:8]
set yrange[-29:-20]
p \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "L_tot", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($39/$13**2)) ti "L_H2", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($34/$13**2)) ti "L_[Cii]158", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($37/$13**2)) ti "L_[Oi ]063", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($38/$13**2)) ti "L_DustRec", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10(($31+$32+$33)/$13**2)) ti "L_rr", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($29/$13**2)) ti "L_res", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($30/$13**2)) ti "L_ff", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($41/$13**2)) ti "H_PI", \
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($42/$13**2)) ti "H_PE",\
"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($45/$13**2)) ti "H_H2"
