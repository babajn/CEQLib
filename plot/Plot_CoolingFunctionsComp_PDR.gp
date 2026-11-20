reset
set multiplot layout 2,2

set key right bottom
set xlabel "log T [K]"
set ylabel "log Cooling/nH^2 [erg cm3/s]"
set xrange[1:8]
set yrange[-32:-20]
p \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+18_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=0.1 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+18_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+18_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=10.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+19_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=0.1 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+19_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+19_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=10.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+20_ISRF1.0e+02_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=0.1 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+20_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=10.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+21_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=0.01 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=0.1 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+21_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+21_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=10.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+02_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=100 Nt=21 G0=0.1 Z=1.0"

p \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+18_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+19_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+21_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=1.0 Z=1.0"

p \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e-01_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=0.1 Nt=20 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+00_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=1.0 Nt=20 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10. Nt=20 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+02_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=100 Nt=20 G0=1.0 Z=1.0"

p \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+18_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=0.01 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+19_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=0.1 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=1.0 Z=1.0", \
"out/PDRCoolingFunction_nH1.0e+01_Nt1.0e+21_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=10. Z=1.0"

#"out/CIECoolingFunction_H2_nH1.0e+01_Nt1.0e+18_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+18_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=0.1 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+18_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=1.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+18_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=18 G0=10.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+19_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=0.1 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+19_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=1.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+19_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=19 G0=10.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+20_ISRF1.0e+02_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=0.1 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+20_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=1.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+20_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=20 G0=10.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-02_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=0.01 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=0.1 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e+00_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=1.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+01_Nt1.0e+21_ISRF1.0e+01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=10 Nt=21 G0=10.0 Z=1.0", \
#"out/PDRCoolingFunction_H2_nH1.0e+02_Nt1.0e+21_ISRF1.0e-01_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "nH=100 Nt=21 G0=0.1 Z=1.0"

