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
def bwa(read, reference, tolerance=0, bwt_data=None, suffix_array=None):
    
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

## case 1: simple case 
read="GCA"
reference="CGATGCACCGGT"
print(bwa(read, reference, tolerance=0))

## case 2: where i failed last time 
read="GCAG"
reference="CGATGCACCGGTACTGGATCGATCGATCGAGTGCTAGCGTAGCGAGGCATGGATCAGGCAG"
print(bwa(read, reference, tolerance=0))


## case 3: no match 
read="GCA"
reference="GACGATGAACCGGTGCCA"
print(bwa(read, reference, tolerance=0))

## case 4: multiple matches. 
read="GCA"
reference="GCACGATGGCAAACCGGTGCGCACA"
print(bwa(read, reference, tolerance=0))

## case 5: read reads and references from the bowtie 2 github https://github.com/BenLangmead/bowtie2
read_1 = 'TGAATGCGAACTCCGGGACGCTCAGTAATGTGACGATAGCTGAAAACTGTACGATAAACNGTACGCTGAGGGCAGAAAAAATCGTCGGGGACATTNTAAAGGCGGCGAGCGCGGCTTTTCCG'
read_2 = 'NTTNTGATGCGGGCTTGTGGAGTTCAGCCGATCTGACTTATGTCATTACCTATGAAATGTGAGGACGCTATGCCTGTACCAAATCCTACAATGCCGGTGAAAGGTGCCGGGATCACCCTGTGGGTTTATAAGGGGATCGGTGACCCCTACGCGAATCCGCTTTCAGACGTTGACTGGTCGCGTCTGGCAAAAGTTAAAGACCTGACGCCCGGCGAACTGACCGCTGAGNCCTATGACGACAGCTATCTCGATGATGAAGATGCAGACTGGACTGC'
read_3 = 'ATCGCCCGCAGACACCTTCACGCTGGACTGTTTCGGCTTTTACAGCGTCGCTTCATAATCCTTTTTCGCCGCCGCCATCAGCGTGTTGTAATCCGCCTGCAGGATTTTCCCGTCTTTCNGTGCCTTGNTCAGTTCTTCCTGACGGGCGGTATATTTCTGCAGCGGCGTCTGCAGCCGTTCGTNAGCCTTCTGCGCCTCTTCGGTATATTTCAGCCGTGACGCTTCGGTATCGCTCCGCTGCTGCGCATTTTTCTCCTCTTGAGTCTGCTGCTCAGCCTTCTTTCGGGCGGCTTCAAGCGCAAGACGGGCCTTTTCACGATCATCCCAGTAACGCGCCC'
ref_file= open(r'/home/buzgalbraith/school/spring_2023/DNA-Mapping-Algorithms/data/lambda_virus.fa', "r")
ref = ref_file.readlines()[1:]
reference = "".join(ref).replace('\n', '')
print(bwa(read_1, reference, tolerance=0))