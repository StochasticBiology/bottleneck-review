set xlabel "b"
set ylabel "V'(h)"
set key bottom left
set logscale y
set yrange [0.001:1]
set label 1 at 0.8, 0.01 "Lower g"
set label 2 at 0.45, 0.8 "Higher g"

plot "simulate-wright.txt" u ($1 == 1 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#FF8800" title "N = 1", "" u ($1 == 2 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#EE8811" title "N = 2", "" u ($1 == 4 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#DD8822" title "N = 4", "" u ($1 == 8 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#CC8833" title "N = 8", "" u ($1 == 16 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#BB8844" title "N = 16", "" u ($1 == 32 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#AA8855" title "N = 32", "" u ($1 == 64 ? $4/0.25 : 1/0):(exp(-($2+1)/$1)):1 w l lw 2 lc rgbcolor "#998866" title "N = 64"    
