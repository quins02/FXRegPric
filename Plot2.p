set term x11
set autoscale
set xlabel "Time (Month)"
set ylabel "Option Value ($)"
set title "Risk Profile of Option over Time\n**From Regressor Pricer**"
plot "tmp/CVA.dat", "tmp/CVAMean.dat" with lines, "tmp/CVAPMean.dat" with lines
# 