set term x11
set autoscale
set xlabel "Rate1"
set ylabel "Rate2"
set zlabel "Option Value ($)"
set title "Interest rate Vs Option Value at Different Time to Expiry\n**From Regressor Pricer**"
splot "tmp/RateV1.dat","tmp/RateV2.dat","tmp/RateV3.dat","tmp/RateV4.dat","tmp/RateV5.dat","tmp/RateV6.dat","tmp/RateV7.dat","tmp/RateV8.dat","tmp/RateV9.dat","tmp/RateV10.dat","tmp/RateV11.dat","tmp/RateV12.dat" 

# 