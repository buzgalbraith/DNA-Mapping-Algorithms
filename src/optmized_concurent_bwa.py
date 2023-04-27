import numpy as np 
from collections import Counter

def generate_suffix_array(reference_array):
    
    suffix_array = np.array(list(sorted(range(reference_array.shape[0]), key=lambda i:reference[i:])))
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
        suffix_array = generate_suffix_array(reference_array)
    bwt_ref = bwt(reference_array, suffix_array)

    rank_map = lf_mapping(bwt_ref,  letters | set([eos]) )

    for i, j in rank_map.items():
        np.append(j, [0])

    return letters, bwt_ref, rank_map, counts, suffix_array

def bwa(read, reference, tolerance=0, bwt_data=None, suffix_array=None):

    results = []
     
    if len(read) == 0:
        return("Empty read")
    if bwt_data is None:
        bwt_data = generate_all(reference, suffix_array=suffix_array)

    letters, bwt, rank_map, count, suffix_array = bwt_data

    if len(letters) == 0:
        return("Empty reference")

    if not set(read).issubset(letters):
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
# case 1: simple case 
read="GCA"
reference="CGATGCACCGGT"
print(bwa(read, reference, tolerance=0))

# case 2: where i failed last time 
read="GCAG"
reference="CGATGCACCGGTACTGGATCGATCGATCGAGTGCTAGCGTAGCGAGGCATGGATCAGGCAG"
print(bwa(read, reference, tolerance=0))


## case 3: no match 
read="GCATTWO"
reference="GCATTWEDTI"
print(bwa(read, reference, tolerance=0))

# ## case 4: multiple matches. 
read="GCA"
reference="GCACGATGGCAAACCGGTGCGCACA"
print(bwa(read, reference, tolerance=0))

# # #case 5: reads and references from the bowtie 2 github https://github.com/BenLangmead/bowtie2
read_1 = 'TGAATGCGAACTCCGGGACGCTCAGTAATGTGACGATAGCTGAAAACTGTACGATAAACNGTACGCTGAGGGCAGAAAAAATCGTCGGGGACATTNTAAAGGCGGCGAGCGCGGCTTTTCCG' ## 10m 5.3 second s
read_2 = 'CGCCAAAAGTGAGAGGCACCTGTCAGATTGAGCGTGCAGCCAGTGAATCCCCGCATTTTATGCGTTTTCATGTTGCCTGCCCGCATTGCGGGGA'

# ref_file= open(r'/home/buzgalbraith/work/school/spring_2023/DNA-Mapping-Algorithms/data/lambda_virus.fa', "r")
# ref = ref_file.readlines()[1:]
# reference = "".join(ref).replace('\n', '')

print(bwa(read_2, reference, tolerance=0)) ## took 7 ish min returns the correct result that there is no match.
