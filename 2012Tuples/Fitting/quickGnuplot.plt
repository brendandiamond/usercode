reset
set nokey
set term aqua enh

F(x)= A*x + B
set ytics nomirror

fit F(x) "BestFitParameters.txt" using 1:2:3 via A,B
set title "MonteCarlo invariant mass from fit of Crystal Ball + P(2)"
set xlabel "Generated m_{{/Symbol g}{/Symbol g}} (GeV)"
set ylabel "Reconstructed m_{{/Symbol g}{/Symbol g}} (GeV)"
set xrange [180:450]
set yrange [120:450]
set y2range [-0.2:1.0]
set y2tics border
plot "BestFitParameters.txt" using 1:2:3  with yerrorbars, F(x), "BestFitParameters.txt" using 1:(F($1) - $2):3 axes x1y2 title "" with yerrorbars

# fit F(x) "BestFitParameters.txt" using 1:4:5 via A,B
# set title "MonteCarlo {/Symbol s} from fit of Crystal Ball + P(2)"
# set xlabel "Generated m_{{/Symbol g}{/Symbol g}} (GeV)"
# set ylabel "{/Symbol s} (GeV)"
# set xrange [180:450]
# set yrange [0:5]
# set y2range [-0.5:2.0]
# set y2tics border
# plot "BestFitParameters.txt" using 1:4:5  with yerrorbars, F(x), "BestFitParameters.txt" using 1:(F($1) - $4):5 axes x1y2 title "" with yerrorbars

XVAL=0.2
YVAL=0.6
set label sprintf("y = %.2E x + %.2E", A, B) at screen XVAL,screen YVAL

set x2zeroaxis lt -1
set y2label "Residuals"
replot

set term post enh color eps
set output "mass.eps"
replot
