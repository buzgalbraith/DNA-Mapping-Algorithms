#!/bin/bash
x=4
while [ $x -le $# ]
do
    if test -f bwa_time_${!x}.txt; then
        rm -r bwa_time_${!x}.txt
    fi

    if test -f bwa_out.txt; then
        rm -r bwa_out.txt
    fi


    repeates=0
    while (( "$repeates" <  "$2"))
    do
        runs=0
   while (( "$runs" <  "$3"))
        do
       # echo  "${!x}"
        { time mpirun  -n $1  --oversubscribe  python3  optimized_paraellel_bwa.py >> bwa_out.txt  ${!x} ;} 2>> bwa_time_${!x}.txt
        runs=$(( $runs + 1 ))
        done 
    repeates=$(( $repeates + 1 ))
    
    done 
x=$(( $x + 1 ))
done
# argumentsnumber of workers,  number of runs,  number of repeates,  then a varible number of arguments for the runs you want to try 