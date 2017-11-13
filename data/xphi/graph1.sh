#!/usr/bin/gnuplot

set term png
set output "h-relations_levels.png"

set xrange [0:2048]
#plot "h2plot.dat" using 1:2 title "Level 2", "h3plot.dat" using 1:2 title "Level 3" , "h4plot.dat" using 1:2 title "Level 4"
#plot "h1plot.dat" using 1:2 title "Level 1", "h2plot.dat" using 1:2 title "Level 2"
plot "h1.dat" using 1:2 title "Level 1" with lines, \
     "h2.dat" using 1:2 title "Level 2" with lines, \
     "h3.dat" using 1:2 title "Level 3" with lines, \
     "h4.dat" using 1:2 title "Level 4" with lines
