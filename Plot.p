set term x11
set autoscale
set xlabel "Time"
set ylabel "FX Rate"
set zlabel "Option Value ($)"
set title "Value Functions for Options at Different Time to Expiry\n**From Regressor Pricer**"
splot "tmp/BucketFit1.dat","tmp/BucketFit2.dat","tmp/BucketFit3.dat","tmp/BucketFit4.dat","tmp/BucketFit5.dat","tmp/BucketFit6.dat","tmp/BucketFit7.dat","tmp/BucketFit8.dat","tmp/BucketFit9.dat","tmp/BucketFit10.dat","tmp/BucketFit11.dat","tmp/BucketFit12.dat" 

# 