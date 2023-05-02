import random 
import os
import re 
import timeit
import numpy as np
import sys
def generate_suffix_array(reference):
    """
    Generates the suffix array of a reference string.

    Parameters:
    -----------
    reference: str
        The reference string to generate the suffix array for.

    Returns:
    --------
    List : A list of indices that represents the suffix array of the input reference string.
    """
    suffix_array = list(sorted(range(len(reference)), key=lambda i:reference[i:]))
    return suffix_array

def bwt(reference, suffix_array=None):
    """
    Generates the Burrows-Wheeler Transform (BWT) of a reference string using its suffix array.

    Parameters:
    -----------
    reference: str
        reference string.
    suffix_array: list
        The suffix array of the reference string.

    Returns:
    --------
    str: the BWT of the input reference string.
    """
    if suffix_array is None:
        suffix_array = generate_suffix_array(reference)
    bwt_ref = "".join(reference[idx - 1] for idx in suffix_array)
    return bwt_ref
def get_all_ranks(bwt_ref, letters=None):
    """
    Returns a dictionary containing the rank of each letter in the Burrows-Wheeler Transform of the reference string.

    Args:
    - bwt_ref (str): The Burrows-Wheeler Transform of the reference string.
    - letters (set, optional): The set of letters to get the ranks for. If not specified, the set of all unique characters in bwt_ref is used.

    Returns:
    - dict: A dictionary where the keys are the letters in letters and the values are lists of integers representing the rank of each occurrence of the letter in bwt_ref.
    """
    if letters is None:
        letters = set(bwt_ref)
        
    result = {letter:[0] for letter in letters}
    result[bwt_ref[0]] = [1]
    for letter in bwt_ref[1:]:
        for i, j in result.items():
            j.append(j[-1] + (i == letter))
    return(result)

def rank(bwt_ref,char,index):
    """
    counts accourances of one charecter up to a certain point in the string. 

    Args:
    - bwt_ref (str): The Burrows-Wheeler Transform of the reference string.
    - char (str): The character to count the occurrences of.
    - index (int): The number of characters to consider from the beginning of bwt_ref.

    Returns:
    - int: The number of occurrences of char in the first index characters of bwt_ref.
    """
    j=0
    occurrences=0
    while(j<=index):
        if(bwt_ref[j]==char):
            occurrences=occurrences+1
        j+=1
    return occurrences 

def char_rank(bwt_ref, char):
    """
    get the rank for for a single charicter at every position in an array

    Args:
    - bwt_ref (str): The Burrows-Wheeler Transform of the reference string.
    - char (str): The character to count the occurrences of.
    Returns:
    - dict: The rank mapping of that charicter
    """
    ranks=[]
    for i in range(len(bwt_ref)):
        ranks.append(rank(bwt_ref, char, i))
    return ranks

def rank_mapping(bwt_ref, letters = None):
    """
    Returns a dictionary containing the rank of each letter in each prefix of bwt_ref.

    Args:
    - bwt_ref (str): The Burrows-Wheeler Transform of the reference string.
    - letters (set, optional): The set of letters to get the ranks for. If not specified, the set of all unique characters in bwt_ref is used.

    Returns:
    - dict: A dictionary where the keys are the letters in letters and the values are lists of integers representing the rank of each occurrence of the letter in each prefix of bwt_ref.
    """
    if letters is None:
        letters = set(bwt_ref)

    fixed_char_rank = lambda char: char_rank(bwt_ref, char)

    result = {letter:fixed_char_rank(letter) for letter in letters}
    return(result)

def get_car_counts(reference):
    """
    counts total occourances of each letter in a reference
    
    Args:
    str: refernce, refernce string
    
    Returns:
    dict: dictionary mapping each letter to its total number of occourances. 
    """
    d={}
    for char in reference:
        try:
            d[char]=d[char]+1
        except:
            d[char]=1
    return d 



def count_occurrences(reference, letters=None):
    """
    Counts the number of occurrences of each character in a reference string.

    Parameters:
    -----------
    reference: str: The reference string.
    letters: set, optional The set of unique characters in the reference string. If None, this will be inferred from the input.

    Returns:
    --------
    dict
        A dictionary mapping each character in the input set to the total number of occurrences in the reference string.
    """
    count = 0
    result = {}
    
    if letters is None:
        letters = set(reference)
        
    c = get_car_counts(reference)
    
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
    suffix_array: list optional
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

    counts = count_occurrences(reference, letters)

    reference = "".join([reference, eos])
    if suffix_array is None:
        suffix_array = generate_suffix_array(reference)
    bwt_ref = bwt(reference, suffix_array)
    rank_map = rank_mapping(bwt_ref, letters | set([eos]))

    for i, j in rank_map.items():
        j.extend([j[-1], 0])

    return letters, bwt_ref, rank_map, counts, suffix_array
def bwa_run(read, reference, bwt_data=None, suffix_array=None):
    """
    Runs Burrows-Wheeler Aligner (BWA) algorithm on a single read and a reference sequence.

    Parameters:
    read (str): A string representing the read to be aligned.
    reference (str): A string representing the reference genome sequence.
    bwt_data (tuple, optional): A tuple containing Burrows-Wheeler Transform (BWT) data 
                                generated from reference genome sequence. Defaults to None.
    suffix_array) lisr, optional: An array containing the suffix array generated
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

    if not set(read) <= letters:
        return []

    length = len(bwt)

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
## testing 

def bwa(ref_path, read_path, num_reads=100):
    """
    bwa function takes a reference file path, a read file path and an optional parameter of number of reads to process. It uses the bwa algorithm to align each read in the read file to the reference sequence and returns a dictionary with the matches for each read.
    
    Args:
    
    ref_path (str): path to the reference file.
    read_path (str): path to the read file.
    num_reads (int): number of reads to process. Default is 100.
    Returns:
    
    out_dict (dict): A dictionary where each key is a read sequence and the corresponding value is a list of matches in the reference sequence. If out is set to False, then a dictionary is returned containing only the matches for each sequence.
    """
    reads = get_n_reads(read_path, num_reads=num_reads)
    ref = get_ref(ref_path)
    out_dict = {}
    for read in reads: 
        out_dict[read] = bwa_run(read, ref)
    return out_dict



def get_n_reads(read_path, num_reads=100):
    """
    Reads a fastq file and extracts a specified number of reads randomly.
    
    Args:
        read_path (str): The path to the fastq file to read.
        num_reads (int, optional): The number of reads to extract. Default is 100.
        
    Returns:
        list: A list of reads with length num_reads, each read being a string of nucleotide bases.
    """
    read_file = open(read_path)
    total_reads = sum(1 for line in read_file) ## count total number of lines 
    reads = []
    i=0
    while i<num_reads:
        read_pos =  random.randrange(total_reads)
        read_file.seek(read_pos)
        read = read_file.readline()
        try:
            read = re.search('^[ACTG]{7,}[ACTG]$', read).group(0)
            reads.append(read)
            i+=1
        except:
            pass
    return reads
def get_ref(ref_path):
    """
    Reads a reference file and extracts the nucleotide bases.
    
    Args:
        ref_path (str): The path to the reference file to read.
        
    Returns:
        str: A string containing only the nucleotide bases of the reference genome.
    """    
    ref_file= open(ref_path, "r")
    ref = ref_file.readlines()
    ref = "".join(re.findall('[ACTG]',"".join(ref).replace('\n', '')))
    return ref
try: ## makes sure that current path is in the data file from our github repo
    os.chdir('..')
    os.chdir('..')
    os.chdir('..')
    os.getcwd()
    os.chdir('./data')  
except: 
    pass 
read_path = 'SRR11528307_R2.fastq'
ref_path = 'SRR11528307_SarsCov2.fna'

num_repeats = int(sys.argv[1])
num_runs = int(sys.argv[2])
read_numbers = list(map(int, sys.argv[3:]))

for num_read in read_numbers:
    timing = timeit.repeat('bwa(ref_path,read_path, num_reads={0})'.format(num_read),globals=globals(),repeat=num_repeats, number=num_runs)   
    timing = list(map(lambda x: x/num_runs, timing)) ## the repeat function returns the sum of the run times so we take the mean of each repeat
    print("for {0} reads\nrun times were {1}\nmean={2}, var={3} ".format(num_read, timing, np.mean(timing), np.var(timing)))
    print("-"*50)
