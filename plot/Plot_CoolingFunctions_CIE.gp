reset
set key right bottom
set xlabel "log T [K]"
set ylabel "log Cooling/nH^2 [erg cm3/s]"
set xrange[1:8]
set yrange[-29:-20]
p \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e-02.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=0.01", \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs3.0e-02.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=0.03", \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e-01.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=0.10", \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs1.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=1.00", \
"out/CIECoolingFunction_nH1.0e+01_Nt1.0e+18_Zs3.0e+00.dat" u (log10($9)):(log10($46/$13**2)) ti "CIE Nt=18 Z=3.00"

