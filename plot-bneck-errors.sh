reset
set xrange [0:50]
set xtics (2, 5, 10, 20, 30, 40, 50)
set yrange [0:100]
set label 1 at 5,95 "..."
unset key
set xlabel "Number of samples"
set ylabel "Estimated bottleneck size N"
plot "bneck-errors.txt" u 1:2:3:4 w yerrorbars lw 2 lc rgbcolor "#FF8888"
