reset
set multiplot
set size 0.33,0.16

set boxwidth 0.05
set style fill solid
set yrange [0.0:*]
unset ytics

set xrange [-0.05 : 1.05]
unset key
set origin 0.33,0.83
plot "simulate-kimura.txt" u 1:2 w boxes
set origin 0.33,0.66
plot "simulate-kimura.txt" u 1:3 w boxes
set origin 0.33,0.5
plot "simulate-kimura.txt" u 1:4 w boxes
set origin 0.33,0.33
plot "simulate-kimura.txt" u 1:5 w boxes
set origin 0.33,0.16
plot "simulate-kimura.txt" u 1:7 w boxes
set origin 0.33,0.0
plot "simulate-kimura.txt" u 1:9 w boxes


set origin 0.0,0.83
plot "simulate-binomial.txt" u 1:2 w boxes
set origin 0.0,0.66
plot "simulate-binomial.txt" u 1:3 w boxes
set origin 0.0,0.5
plot "simulate-binomial.txt" u 1:4 w boxes
set origin 0.0,0.33
plot "simulate-binomial.txt" u 1:5 w boxes
set origin 0.0,0.16
plot "simulate-binomial.txt" u 1:7 w boxes
set origin 0.0,0.0
plot "simulate-binomial.txt" u 1:9 w boxes

set origin 0.66,0.83
# beta distribution model not well defined for N = 1
set origin 0.66,0.66
plot "simulate-beta.txt" u 1:3 w boxes
set origin 0.66,0.5
plot "simulate-beta.txt" u 1:4 w boxes
set origin 0.66,0.33
plot "simulate-beta.txt" u 1:5 w boxes
set origin 0.66,0.16
plot "simulate-beta.txt" u 1:7 w boxes
set origin 0.66,0.0
plot "simulate-beta.txt" u 1:9 w boxes

