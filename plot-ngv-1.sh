reset
set multiplot
set size 0.8,1
set origin 0.2,0
set logscale x
set logscale x2
set logscale y
set key bottom left
set xrange [1:500]
set xlabel "N"
set ylabel "V'(h)"

set xtics nomirror
set xtics (1, 2, 5, 10, 20, 50, 100, 200, 500)

# this simply plots (1 - (1-1/x)^g) for different g
plot (1.-(1-1./x)**1) lw 2 lc rgbcolor "#FF8800" title "g = 1", (1.-(1-1./x)**2) lw 2 lc rgbcolor "#EE8811" title "g = 2", (1.-(1-1./x)**5) lw 2 lc rgbcolor "#DD8822" title "g = 5", (1.-(1-1./x)**10) lw 2 lc rgbcolor "#CC8833" title "g = 10", (1.-(1-1./x)**20) lw 2 lc rgbcolor "#BB8844" title "g = 20", (1.-(1-1./x)**50) lw 2 lc rgbcolor "#AA8855" title "g = 50", (1.-(1-1./x)**100) lw 2 lc rgbcolor "#998866" title "g = 100", (1.-(1-1./x)**200) lw 2 lc rgbcolor "#888877" title "g = 200", (1.-(1-1./x)**500) lw 2 lc rgbcolor "#778888" title "g = 500", (1.-(1-1./x)**1000) lw 2 lc rgbcolor "#668899" title "g = 1000"

set size 0.2
set border 1
unset logscale
set xrange [0:1]
unset xtics
unset x2tics
unset ytics
unset xlabel
unset ylabel
unset x2label
unset key
set xtics (0,1) nomirror

set style fill solid

# kimura distributions (and final normal approximation) for illustration of corresponding distribution widths

set origin 0,0.8
plot "simulate-kimura.txt" u 1:2 w boxes

set origin 0,0.55
plot "simulate-kimura.txt" u 1:5 w boxes

set origin 0,0.3
plot "simulate-kimura.txt" u 1:8 w boxes

set origin 0,0.05
sigma = sqrt(0.001*0.25)
plot 1./sqrt(2*3.1416*sigma**2)*exp(-(x-0.5)**2/(2.*sigma**2)) w filledcu
