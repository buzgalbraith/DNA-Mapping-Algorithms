#!/bin/bash

# argumentsnumber of workers,  number of runs,  number of repeates,  then a varible number of arguments for the runs you want to try 

x=4 # optinal argumnets (ie number of reads we want to try start at comand line arg 4)
while [ $x -le $# ] ## while x is less than or equal to the number of argumnets left
do
    if test -f bwa_time_${!x}.txt; then ## if a file with the name we are planning on outputting our times too exists delete it
        rm -r bwa_time_${!x}.txt
    fi

    if test -f bwa_out.txt; then  ## if a file with the name we are planning on outputting our dictionary too exists delete it
        rm -r bwa_out.txt
    fi

    ## does a while loop for each repeate and run arugments taken from comand line 
    repeates=0 
    while (( "$repeates" <  "$2"))
    do
        runs=0 
   while (( "$runs" <  "$3"))
        do
       # echo  "${!x}"
        { time mpirun  -n $1  --oversubscribe  python3  optimized_paraellel_bwa.py >> bwa_out.txt  ${!x} ;} 2>> bwa_time_${!x}.txt
        ## passes the arguments and runs optimized_paraellel_bwa.py. append the resulting dict to the bwa_out file
        ## append the timing information to a file named bwa_time_run_number.txt 
        runs=$(( $runs + 1 ))
        done 
    repeates=$(( $repeates + 1 ))
    
    done 
x=$(( $x + 1 ))
done
