#!/bin/bash
# Script for running different tests

clear

if [ $# -eq 0 ]; then # Run the tests
	 
	echo "\nTests starting...\n"
	rm -r "./data"

	num_tests=$(ls bin/ | wc -l)
	nr_iterations=10

	for iter in $(seq 1 $nr_iterations); do
		echo "======= Iteration $iter / $nr_iterations ======= \n"
		current_test=1 
		mkdir -p "./data/$iter"
		for file in bin/* ; do
			echo "======= $current_test / $num_tests : $file ======= \n"
			cp $file "conv_check" 
			
			./conv_check | tee "./data/$iter/$(basename $file).dat"

			current_test=$(($current_test+1))
		done
	done
else # Compile the binarys
echo "Compiling binarys"
	num_tests=$(ls tests/*.f90 | wc -l)
	current_test=1 
	mkdir -p "./bin/"
	for file in tests/*.f90 ; do
		echo "======= $current_test / $num_tests : $file ======= \n"
		cp $file "papi_counters.f90" 
		touch "conv_check.f90"

		make conv_check
		
		mv "conv_check" "./bin/$(basename $file .f90)"

		current_test=$(($current_test+1))
	done
echo "Compiling Done!"
fi

