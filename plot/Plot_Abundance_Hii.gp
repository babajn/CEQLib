reset
set key left bottom
set xlabel "log T [K]"
set ylabel "log Fraction"
set xrange[2:6]
set yrange[-5:0]
p \
"out/PostProcessMap.dat" u (log10($9)):(log10($17/$13)) ti "e/H", \
"out/PostProcessMap.dat" u (log10($9)):(log10($18/$13)) ti "HII/H", \
"out/PostProcessMap.dat" u (log10($9)):(log10($19/$13)) ti "HI/H", \
"out/PostProcessMap.dat" u (log10($9)):(log10($20/$14)) ti "HeIII/He", \
"out/PostProcessMap.dat" u (log10($9)):(log10($21/$14)) ti "HeII/He", \
"out/PostProcessMap.dat" u (log10($9)):(log10($22/$14)) ti "HeI/He", \
"out/PostProcessMap.dat" u (log10($9)):(log10($23/$15)) ti "CIII/C", \
"out/PostProcessMap.dat" u (log10($9)):(log10($24/$15)) ti "CII/C", \
"out/PostProcessMap.dat" u (log10($9)):(log10($25/$15)) ti "CI/C", \
"out/PostProcessMap.dat" u (log10($9)):(log10($26/$16)) ti "OIII/O", \
"out/PostProcessMap.dat" u (log10($9)):(log10($27/$16)) ti "OII/O", \
"out/PostProcessMap.dat" u (log10($9)):(log10($28/$16)) ti "OI/O"
