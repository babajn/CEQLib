reset
set key left
set xlabel "log n [H/cc]"
set ylabel "log T [K], p/k [K/cc]"
set xrange[-3:+3]
set yrange[0:6]
p "out/PostProcessMap.dat" u (log10($13)):(log10($9)) ti "n-T", "out/PostProcessMap.dat" u (log10($13)):(log10($9*$13)) ti "n-P"
