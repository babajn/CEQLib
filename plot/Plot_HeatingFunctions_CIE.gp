reset
set key left
set xlabel "log T [K]"
set ylabel "log Heating/nH^2 [erg cm3/s]"
set xrange[2:6]
set yrange[-29:-20]
p \
"out/CIECoolingFunction_Nt1.0e+18_ISRF0.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "CIE Nt=18 Z=1.0", \
"out/CoolingFunction_Nt1.0e+18_ISRF0.1_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=18 G0=0.1 Z=1.0", \
"out/CoolingFunction_Nt1.0e+18_ISRF1.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=18 G0=1.0 Z=1.0", \
"out/CoolingFunction_Nt1.0e+18_ISRF10.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=18 G0=10.0 Z=1.0", \
"out/CoolingFunction_Nt1.0e+19_ISRF0.1_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=19 G0=0.1 Z=1.0", \
"out/CoolingFunction_Nt1.0e+19_ISRF1.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=19 G0=1.0 Z=1.0", \
"out/CoolingFunction_Nt1.0e+19_ISRF10.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=19 G0=10.0 Z=1.0", \
"out/CoolingFunction_Nt1.0e+20_ISRF0.1_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=20 G0=0.1 Z=1.0", \
"out/CoolingFunction_Nt1.0e+20_ISRF1.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=20 G0=1.0 Z=1.0", \
"out/CoolingFunction_Nt1.0e+20_ISRF10.0_Zs1.0e+00.dat" u (log10($9)):(log10($47/$13**2)) ti "Nt=20 G0=10.0 Z=1.0"

