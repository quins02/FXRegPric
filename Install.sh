#! /bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # 
hash gsl-config 2>/dev/null || { echo >&2 "This program requires the GSL libraries. Please install them. Aborting"; exit 1; }
echo "Compiling Regression Pricing Model In Place..."
mkdir tmp >/dev/null 2>&1
mkdir bin >/dev/null 2>&1
CF="$(gsl-config --cflags --libs)"
g++ -Wall src/VProfile.cpp src/mat_opp.cpp src/opt_eval.cpp src/rng_func.cpp src/download.cpp -O3 -o FXRegProfile.exe -std=c++11 -lcurl ${CF}
echo "Required DLLs for portability:"
for i in $(cygcheck ./FXRegProfile.exe 2>/dev/null ); 
	do 
		if [[ ! "$i"  =~ "cygcheck" ]] 
			then
				cp $i bin
				echo Copied to .../bin: $i 
		fi 
	done
echo "Success"
echo "Compiled as Regressor.exe."