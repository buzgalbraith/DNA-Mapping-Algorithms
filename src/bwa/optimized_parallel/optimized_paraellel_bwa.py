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
num_reads = int(sys.argv[1]) ## takes the number of reads as a command line arg
comm = MPI.COMM_WORLD 
num_workers = comm.Get_size() ## save the size as the number of workers
rank = comm.Get_rank()


def generate_suffix_array(reference):
    """
    Generates the suffix array of a reference string.

    Parameters:
    -----------
    reference: str
        The reference string to generate the suffix array for.

    Returns:
    --------
    numpy.ndarray : An array of indices that represents the suffix array of the input reference string.
    """
    suffix_array = np.array(list(sorted(range(len(reference)), key=lambda i:reference[i:])))
    return suffix_array.astype("int")

def bwt(reference_array, suffix_array):
    """
    Generates the Burrows-Wheeler Transform (BWT) of a reference string using its suffix array.

    Parameters:
    -----------
    reference_array: numpy.ndarray
        The array of characters in the reference string.
    suffix_array: numpy.ndarray
        The suffix array of the reference string.

    Returns:
    --------
    numpy.ndarray the BWT of the input reference string.
    """
    bwt_ref_array= reference_array[suffix_array-1]
    return bwt_ref_array

def lf_mapping(bwt_ref, letters):
    """
    Generates the last-first mapping (LF mapping) of the BWT of a reference string.

    Parameters:
    -----------
    bwt_ref: numpy.ndarray
        The BWT array  of the reference string.
    letters: set
        The set of unique characters in the reference string.

    Returns:
    --------
    dict
        A dictionary mapping each character in the input set to an array of indices in the BWT.
    """
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
    """
    Counts the number of occurrences of each character in a reference string.

    Parameters:
    -----------
    reference_array: numpy.ndarray
        The array of characters in the reference string.
    letters: set, optional
        The set of unique characters in the reference string. If None, this will be inferred from the input.

    Returns:
    --------
    dict
        A dictionary mapping each character in the input set to the total number of occurrences in the reference string.
    """
    count = 0
    result = {}
    c = Counter(reference_array)

    for letter in sorted(letters):
        result[letter] = count
        count += c[letter]

    return result
def update(begin, end, char, rank_map, counts):
    """
    Updates the beginning and ending indices for a partial match based on a single character.

    Parameters:
    -----------
    begin: int
        The beginning index of the partial match.
    end: int
        The ending index of the partial match.
    char: str
        The character to use for the update.
    rank_map: dict
        A dictionary mapping each character to an array of indices in the BWT.
    counts: dict
        A dictionary mapping each character to the total number of occurrences in the reference string.

    Returns:
    --------
    tuple
        A tuple of two integers representing the new beginning and ending indices for the partial match.
    """
    beginning = counts[char] + rank_map[char][begin - 1] + 1
    ending = counts[char] + rank_map[char][end]

    return(beginning,ending)

def generate_all(reference, suffix_array=None, eos="$"):
    """
    Generates all data structures needed for BWA alignment.

    Parameters:
    -----------
    reference: str
        The reference string to generate the data structures for.
    suffix_array: numpy.ndarray, optional
        The suffix array of the reference string. If None, this will be generated automatically.
    eos: str, optional
        The end-of-string character to add to the end of the reference string.

    Returns:
    --------
    tuple
        A tuple of five items:
            1. A set of unique characters in the reference string.
            2. The BWT of the reference string.
            3. A dictionary mapping each character to an array of indices in the BWT.
            4. A dictionary mapping each character to the total number of occurrences in the reference string.
            5. The suffix array of the reference string.
    """
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

def bwa_run(read, reference, bwt_data=None, suffix_array=None):
    """
    Runs Burrows-Wheeler Aligner (BWA) algorithm on a single read and a reference sequence.

    Parameters:
    read (str): A string representing the read to be aligned.
    reference (str): A string representing the reference genome sequence.
    bwt_data (tuple, optional): A tuple containing Burrows-Wheeler Transform (BWT) data 
                                generated from reference genome sequence. Defaults to None.
    suffix_array (numpy.ndarray, optional): A NumPy array containing the suffix array generated
                                            from reference genome sequence. Defaults to None.

    Returns:
    list: A list of integers representing the positions of read in the reference genome sequence.
          Returns an empty list if the read contains characters that are not in the reference genome 
          sequence, or if either the read or the reference genome sequence is empty.

    """
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

    class match_dict(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    working_dict = [match_dict(read=read, begin=0, end=len(bwt)-1)]


    i=0
    while len(working_dict) > 0:
        i=i+1
        p = working_dict.pop()


        read = p.read[:-1]
        last = p.read[-1]

        all_letters = [last] 

        for letter in all_letters:
            begin, end = update(p.begin, p.end, letter, rank_map, count)

            if begin <= end:
                if len(read) == 0:

                    results.extend(suffix_array[begin : end + 1])
                
                    working_dict.append(match_dict(read=read, begin=begin,
                                            end=end))
    return sorted(set(results))
def get_n_reads(read_path, ref, num_reads=100, num_workers=10):
    """
    generates a list of dictionaries where the number of reads is evenly divided across all workers
    and each dictionary has a set of reads as keys with empty list as values, and a key ref with the
    reference string as value. This generates the buffer which is scattered by mpi4py.
    
    Parameters:
    read_path (str): Path to the read 
    num_reads (int, optional): An integer representing the number of reads to be aligned.
                               Defaults to 100.
    num_workers (int, optional): An integer representing the number of workers used by mpi4py
                               Defaults to 10.
    Returns:
    dict: a list of dictionaries where the number of reads is evenly divided across all workers
            and each dictionary has a set of reads as keys with empty list as values, and a key ref with the
            reference string as value.
  
    """
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
        match = re.search('^[ACTG]{7,}[ACTG]$', read) 
        if(match!=None):
            l_read = match.group(0)
            chunk_read_list.append({l_read:[]})
            i+=1
            j=j+1
            if(j>=read_per_worker):
                reads.append(chunk_read_list)
                j=0
                chunk_read_list = []  
    return reads
def get_ref(ref_path):
    """
    reads the reference string from file
    
    Parameters:
    ref_path (str): file path to reference genome ,

    Returns:
    str : A string containing the reference genome.
    
    """


    ref_file= open(ref_path, "r")
    ref = ref_file.readlines()
    ref = "".join(re.findall('[ACTG]',"".join(ref).replace('\n', '')))
    return ref
def bwa(ref, read_dict, num_reads=100):
    """
    Runs the BWA algorithm on multiple reads and a reference genome sequence.

    Parameters:
    reference (str): A string representing the reference genome sequence..
    read_dict (dict): An array with all a set of reads as keys and empty list as values. This is the buffer for mpi4py
    num_reads (int, optional): An integer representing the number of reads to be aligned.
                               Defaults to 100.

    Returns:
    dict: A dictionary with each key representing a read and the corresponding value
          representing the positions of the read in the reference genome sequence.

    """
    out_dict = {}
    read_dict = dict(pair for d in read_dict for pair in d.items())
    for read_key in read_dict.keys():
        out_dict[read_key] = bwa_run(read_key, ref)
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
    read_dict_data = get_n_reads(read_path, num_reads=num_reads, ref=ref, num_workers=num_workers) ## initialize buffer to scatter
    
else:
    read_dict_data =  None ## initialize buffer to receive scatter

read_dict_data_chunk = comm.scatter(read_dict_data, root=0) ## scatter section of list with read dictionary to each worker.


read_dict_data_chunk = bwa(ref, read_dict=read_dict_data_chunk, num_reads=num_reads) ## each worker runs bwa algorithm on chunk of read_dict
reads_dict = comm.gather(read_dict_data_chunk, root=0) ## gather the read dictionary back to a list
if rank==0:
    reads_dict = dict(pair for d in reads_dict for pair in d.items()) ## put all elements in the list into a single dictionary 
    print(reads_dict) ## print the dictionary (output will be directed to a file)
