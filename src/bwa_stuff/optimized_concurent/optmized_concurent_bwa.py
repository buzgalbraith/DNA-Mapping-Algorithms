import numpy as np 
from collections import Counter
import random
import re 
import os
import timeit


import sys 
def generate_suffix_array(reference):
    
    suffix_array = np.array(list(sorted(range(len(reference)), key=lambda i:reference[i:])))
    return suffix_array.astype("int")

def bwt(reference_array, suffix_array):
    bwt_ref_array= reference_array[suffix_array-1]
    return bwt_ref_array

def lf_mapping(bwt_ref, letters):
    result = {letter: np.zeros(bwt_ref.shape[0] , dtype=np.int32) for letter in letters}
    counts = {letter: np.zeros(bwt_ref.shape[0] , dtype=np.int32) for letter in letters}

    # initialize the first column of the count matrix
    for letter, count in counts.items():
        count[0] = np.sum(bwt_ref == letter)

    # compute the count matrix using NumPy cumulative sum
    for letter in letters:
        mask = bwt_ref == letter
        counts[letter] = np.cumsum(mask)

    # compute the rank map using the count matrix
    for letter in letters:
        result[letter] = np.concatenate((counts[letter], np.array([0])))
    return result





def count_occurrences(reference_array, letters=None, ):
    count = 0
    result = {}
    c = Counter(reference_array)

    for letter in sorted(letters):
        result[letter] = count
        count += c[letter]

    return result
def update(begin, end, char, rank_map, counts):
    beginning = counts[char] + rank_map[char][begin - 1] + 1
    ending = counts[char] + rank_map[char][end]

    return(beginning,ending)

def generate_all(reference, suffix_array=None, eos="$"):

    letters = set(reference)

    assert eos not in letters,  "{0} already in input string".format(eos)
    reference = "".join([reference, eos])
    reference_array = np.array(list(reference))
    counts = count_occurrences(reference_array[:-1], letters)
    if suffix_array is None:
        suffix_array = generate_suffix_array(reference)
    bwt_ref = bwt(reference_array, suffix_array)
    rank_map = lf_mapping(bwt_ref,  letters | set([eos]) )
    for k in rank_map.keys():
         rank_map[k]=np.hstack([rank_map[k],[rank_map[k][-1]],[0]])
    #return rank_map
    return letters, bwt_ref, rank_map, counts, suffix_array

def bwa_run(read, reference, tolerance=0, bwt_data=None, suffix_array=None):
   # print(reference)
    results = []
     
    if len(read) == 0:
        return("Empty read")
    if bwt_data is None:
        bwt_data = generate_all(reference, suffix_array=suffix_array)
    #return bwt_data
    letters, bwt, rank_map, count, suffix_array = bwt_data

    if len(letters) == 0:
        return("Empty reference")

    if not set(read) <= set(letters):
        return []

    length = bwt.shape[0]

    class Fuzzy(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    fuz = [Fuzzy(read=read, begin=0, end=len(bwt)-1, tolerance=tolerance)]


    i=0
    while len(fuz) > 0:
        i=i+1
        p = fuz.pop()


        read = p.read[:-1]
        last = p.read[-1]

        all_letters = [last] if p.tolerance == 0 else letters

        for letter in all_letters:
            begin, end = update(p.begin, p.end, letter, rank_map, count)

            if begin <= end:
                if len(read) == 0:

                    results.extend(suffix_array[begin : end + 1])
                else:
                    dist = p.tolerance
                    if letter != last:
                
                        dist = max(0, p.tolerance - 1)
                    fuz.append(Fuzzy(read=read, begin=begin,
                                            end=end, tolerance=dist))
  #  print(results)
    return sorted(set(results))
def bwa(ref_path, read_path, num_reads=100, tol=0):
    """ 
    out is for if you want, to print out the results, this is more for validating. 
    so if it is set to False, then the algorithm will just construct a dict that returns the matches for each sequence. 
    """
    reads = get_n_reads(read_path, num_reads=num_reads)
    ref = get_ref(ref_path)
    out_dict = {}
    for read in reads: 
        out_dict[read] = bwa_run(read, ref, tolerance=tol)
    return out_dict


def get_n_reads(read_path, num_reads=100):
    read_file = open(read_path)
    total_reads = sum(1 for line in read_file) ## count total number of lines 
    reads = []
    i=0
    while i<num_reads:
        read_pos =  random.randrange(total_reads)
        read_file.seek(read_pos)
        read = read_file.readline()
        try:
            read = re.search('^[ACTG]{3,}[ACTG]$', read).group(0)
            reads.append(read)
            i+=1
        except:
            pass
    return reads
def get_ref(ref_path):
    ref_file= open(ref_path, "r")
    ref = ref_file.readlines()
    ref = "".join(re.findall('[ACTG]',"".join(ref).replace('\n', '')))
    return ref

try: ## makes sure that current path is in the data file from our github repo
    os.chdir('..')
    os.chdir('..')
    os.getcwd()
    os.chdir('./data')  
except: 
    pass 
read_path = 'SRR11528307_R2.fastq'
ref_path = 'SRR11528307_SarsCov2.fna'
#read_positions_dict = bwa(ref_path,read_path, num_reads=40)


num_repeats = int(sys.argv[1])
num_runs = int(sys.argv[2])
read_numbers = list(map(int, sys.argv[3:]))

for num_read in read_numbers:
    timing = timeit.repeat('bwa(ref_path,read_path, num_reads={0})'.format(num_read),globals=globals(),repeat=num_repeats, number=num_runs)   
    timing = list(map(lambda x: x/num_runs, timing)) ## the repeat function returns the sum of the run times so we take the mean of each repeat
    print("for {0} reads\nrun times were {1}\nmean={2}, var={3} ".format(num_read, timing, np.mean(timing), np.var(timing)))
    print("-"*50)