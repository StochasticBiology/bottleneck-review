reset
set multiplot
set size 0.5

set yrange [1:*]
set logscale

set xlabel "True bottleneck size"
set ylabel "Inferred bottleneck size"

noff = 16
set origin 0,0.5
sel = 1
set label 1 at graph 0.05, 0.9 "s = 1"
plot "sampling-scheme.txt" u ($2 == sel && $3 == noff ? $1 : 1/0):8:10:11 w yerrorbars lw 2 lc rgbcolor "#8888FF" title "Variance", "" u ($2 == sel && $3 == noff ? $1 : 1/0):4:6:7 w yerrorbars lw 2 lc rgbcolor "#FF8888" title "MSD", x notitle 

set origin 0.5,0.5
sel = 0.909091
set label 1 at graph 0.05, 0.9 "s = 0.91"
plot "sampling-scheme.txt" u ($2 == sel && $3 == noff ? $1 : 1/0):8:10:11 w yerrorbars lw 2 lc rgbcolor "#8888FF" title "Variance", "" u ($2 == sel && $3 == noff ? $1 : 1/0):4:6:7 w yerrorbars lw 2 lc rgbcolor "#FF8888" title "MSD", x notitle 

set origin 0,0.
sel = 0.751315
set label 1 at graph 0.05, 0.9 "s = 0.75"
plot "sampling-scheme.txt" u ($2 == sel && $3 == noff ? $1 : 1/0):8:10:11 w yerrorbars lw 2 lc rgbcolor "#8888FF" title "Variance", "" u ($2 == sel && $3 == noff ? $1 : 1/0):4:6:7 w yerrorbars lw 2 lc rgbcolor "#FF8888" title "MSD", x notitle 

set origin 0.5,0.
sel = 0.513158
set label 1 at graph 0.05, 0.9 "s = 0.51"
plot "sampling-scheme.txt" u ($2 == sel && $3 == noff ? $1 : 1/0):8:10:11 w yerrorbars lw 2 lc rgbcolor "#8888FF" title "Variance", "" u ($2 == sel && $3 == noff ? $1 : 1/0):4:6:7 w yerrorbars lw 2 lc rgbcolor "#FF8888" title "MSD", x notitle 
