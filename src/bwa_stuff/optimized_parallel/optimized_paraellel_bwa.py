import numpy as np 
from collections import Counter
import random
import re 
import os
import timeit
from mpi4py import MPI
import sys
# print('Argument List:', str(sys.argv))
read_path = 'SRR11528307_R2.fastq'
ref_path = 'SRR11528307_SarsCov2.fna'
num_reads = int(sys.argv[1])
comm = MPI.COMM_WORLD
num_workers = comm.Get_size()
tol = 0
rank = comm.Get_rank()
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
    results = []
     
    if len(read) == 0:
        return("Empty read")
    if bwt_data is None:
        bwt_data = generate_all(reference, suffix_array=suffix_array)
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
    return sorted(set(results))
def get_n_reads(read_path, ref, num_reads=100, num_workers=10):
    read_file = open(read_path)
    total_reads = sum(1 for line in read_file) ## count total number of lines 
    reads = []
    i=0
    chunk_read_list = []
    j=0
    read_per_worker = num_reads//num_workers
    while i<num_reads:
        read_pos =  random.randrange(total_reads)
        read_file.seek(read_pos)
        read = read_file.readline()       
        #print(j)    
        match = re.search('^[ACTG]{3,}[ACTG]$', read) 
        if(match!=None):
            l_read = re.search('^[ACTG]{7,}[ACTG]$', read).group(0)
            chunk_read_list.append({l_read:[]})
            i+=1
            j=j+1
            if(j>=read_per_worker):
                reads.append(chunk_read_list)
                j=0
                chunk_read_list = []  
    return reads
def get_ref(ref_path):
    ref_file= open(ref_path, "r")
    ref = ref_file.readlines()
    ref = "".join(re.findall('[ACTG]',"".join(ref).replace('\n', '')))
    return ref
def bwa(ref, read_dict, num_reads=100, tol=0):
    out_dict = {}
    read_dict = dict(pair for d in read_dict for pair in d.items())
    for read_key in read_dict.keys():
        out_dict[read_key] = bwa_run(read_key, ref, tolerance=tol)
    return out_dict



try: ## makes sure that current path is in the data file from our github repo
    os.chdir('..')
    os.chdir('..')
    os.getcwd()
    os.chdir('./data')  
except: 
    pass 
ref = get_ref(ref_path)
if rank==0:
    read_dict_data = get_n_reads(read_path, num_reads=num_reads, ref=ref, num_workers=num_workers)
    
else:
    read_dict_data =  None 

read_dict_data_chunk = comm.scatter(read_dict_data, root=0)


read_dict_data_chunk = bwa(ref, read_dict=read_dict_data_chunk, num_reads=num_reads, tol=tol)
reads_dict = comm.gather(read_dict_data_chunk, root=0)
if rank==0:
    reads_dict = dict(pair for d in reads_dict for pair in d.items())
    print(reads_dict)
