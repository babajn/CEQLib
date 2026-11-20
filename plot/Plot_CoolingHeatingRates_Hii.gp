reset
set key left
set xlabel "log T [K]"
set ylabel "log Cooling,Heating/nH^2 [erg cm3/s]"
set xrange[2:6]
set yrange[-28:-20]
p \
"out/PostProcessMap.dat" u (log10($9)):(log10($46/$13**2)) ti "L_tot", \
"out/PostProcessMap.dat" u (log10($9)):(log10($47/$13**2)) ti "H_tot", \
"out/PostProcessMap.dat" u (log10($9)):(log10($34/$13**2)) ti "L_[Cii]158", \
"out/PostProcessMap.dat" u (log10($9)):(log10($37/$13**2)) ti "L_[Oi ]063", \
"out/PostProcessMap.dat" u (log10($9)):(log10($38/$13/$17)) ti "L_DustRec", \
"out/PostProcessMap.dat" u (log10($9)):(log10(($31+$32+$33)/$13**2)) ti "L_rr", \
"out/PostProcessMap.dat" u (log10($9)):(log10($29/$13**2)) ti "L_res", \
"out/PostProcessMap.dat" u (log10($9)):(log10($30/$13**2)) ti "L_ff", \
"out/PostProcessMap.dat" u (log10($9)):(log10($41/$13**2)) ti "H_PI", \
"out/PostProcessMap.dat" u (log10($9)):(log10($42/$13**2)) ti "H_PE",\
"out/PostProcessMap.dat" u (log10($9)):(log10($50/$13**2)) ti "H_alpha"

#p \
#"out/PostProcessMap.dat" u (log10($9)):(log10($46/$13**2)) ti "L_tot", \
#"out/PostProcessMap.dat" u (log10($9)):(log10($47/$13**2)) ti "H_tot"
