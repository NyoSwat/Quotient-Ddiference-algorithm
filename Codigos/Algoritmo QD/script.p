set autoscale
unset log
unset label
set xtic auto
set ytic auto
set grid
set title "Polinomio"
set xlabel "x"
set ylabel "f(x)"
plot 0 title "Cero" with lines,\
"Polinomio.dat" using 1:2 title "Polinomio" with lines,\
 "RaicesReales.dat" using 1:2 title "Raiz" with points
 
