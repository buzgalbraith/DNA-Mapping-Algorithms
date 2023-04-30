import random 
import os
import re 
def generate_suffix_array(reference):
    suffix_array = list(sorted(range(len(reference)), key=lambda i:reference[i:]))
    return suffix_array

def bwt(reference, suffix_array=None):
    if suffix_array is None:
        suffix_array = generate_suffix_array(reference)
    bwt_ref = "".join(reference[idx - 1] for idx in suffix_array)
    return bwt_ref
def get_all_ranks(bwt_ref, letters=None):
    if letters is None:
        letters = set(bwt_ref)
        
    result = {letter:[0] for letter in letters}
    result[bwt_ref[0]] = [1]
    for letter in bwt_ref[1:]:
        for i, j in result.items():
            j.append(j[-1] + (i == letter))
    return(result)

def rank(bwt_ref,char,index):
    j=0
    occurrences=0
    while(j<=index):
        if(bwt_ref[j]==char):
            occurrences=occurrences+1
        j+=1
    return occurrences 

def char_rank(bwt_ref, char):
    ranks=[]
    for i in range(len(bwt_ref)):
        ranks.append(rank(bwt_ref, char, i))
    return ranks

def rank_mapping(bwt_ref, letters = None):
    if letters is None:
        letters = set(bwt_ref)

    fixed_char_rank = lambda char: char_rank(bwt_ref, char)

    result = {letter:fixed_char_rank(letter) for letter in letters}
    return(result)

def get_car_counts(reference):
    d={}
    for char in reference:
        try:
            d[char]=d[char]+1
        except:
            d[char]=1
    return d 



def count_occurrences(reference, letters=None):
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
    beginning = counts[char] + rank_map[char][begin - 1] + 1
    ending = counts[char] + rank_map[char][end]
    return(beginning,ending)



def generate_all(reference, suffix_array=None, eos="$"):
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
def bwa_run(read, reference, tolerance=0, bwt_data=None, suffix_array=None):
    
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

    class Fuzzy(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    fuz = [Fuzzy(read=read, begin=0, end=len(bwt)-1, tolerance=tolerance)]

    while len(fuz) > 0:
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
## testing 

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
    os.getcwd()
    os.chdir('./data')  
except: 
    pass 
read_path = 'SRR11528307_R2.fastq'
ref_path = 'SRR11528307_SarsCov2.fna'

bwa(ref_path,read_path, num_reads=1)
