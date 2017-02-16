set terminal epslatex size 3.5, 2.62 color colortext
set output 'en.tex'

set xlabel '$T\, [J/\mathrm{k}_{\scriptscriptstyle{\rm B}}]$'
set ylabel '$\langle E\rangle/ L^2\, [J]$'
unset key
plot 'energyvsbeta.dat' u 1:2:3 w yerrorbars notitle, ''  u 1:2 w p ps .5 pt 7 notitle
