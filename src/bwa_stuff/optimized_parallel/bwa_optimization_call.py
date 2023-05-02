import subprocess
import os 
import re 
import numpy as np
import sys 
## save the command line arguments to variables to make code more readable
number_of_workers = sys.argv[1]
number_of_repeats = sys.argv[2]
number_of_runs = sys.argv[3]
read_numbers_str = sys.argv[4:]
read_numbers_str = " ".join(read_numbers_str) ## variable number of number of reads we want to try all saved as a string. 

bashCommand = "bash bwa_paralell.bash {0} {1} {2} {3}".format(number_of_workers, number_of_repeats, number_of_runs, read_numbers_str) ## construct a bash command 
subprocess.call(bashCommand, shell=True) ## run the bash command on bwa_paralell.bash 


## here we are getting the times and outputting results 
path = os.getcwd()+"/bwa_time_"  ## change path to where the results are saved
read_numbers = list(map(int, sys.argv[4:])) ## make list of read numbers 
for read_number in read_numbers:
    file = open(path+str(read_number)+".txt") ## open corresponding time file 
    times_list = []
    for line in file: ## use regex to find the "real time numbers"
        if "real" in line:
            times = re.findall('\d{1,}', line)
            time =int(times[0])*60+int(times[1])+int(times[2])/1000 ## convert to seconds 
            times_list.append(time)
    print("for {0} reads\nrun times were {1}\nmean={2}, var={3} ".format(read_number, times_list, np.mean(times_list), np.var(times_list)))
    print("-"*50)

