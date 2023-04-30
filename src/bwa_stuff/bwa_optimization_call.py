import subprocess
import os 
import re 
import numpy as np
number_of_workers = 5
number_of_runs = 3
number_of_repeats = 2
read_numbers = [5, 10]

read_numbers_str = list(map(str, read_numbers))
read_numbers_str = " ".join(read_numbers_str)

bashCommand = "bash bwa_paralell.bash {0} {1} {2} {3}".format(number_of_workers, number_of_repeats, number_of_runs, read_numbers_str)
subprocess.call(bashCommand, shell=True)

path = os.getcwd()+"/bwa_time_"
for read_number in read_numbers:
    file = open(path+str(read_number)+".txt")
    times_list = []
    for line in file:
        if "real" in line:
            times = re.findall('\d{1,}', line)
            time =int(times[0])*60+int(times[1])+int(times[2])/1000
            times_list.append(time)
    print("For {0} reads".format(read_number))
    print("Times were {0}".format(times_list))
    print("mean={0}, var={1}".format(np.mean(times_list), np.var(times_list)))
    print("-"*50)
