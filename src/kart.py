class LMEM:
    """
    An LMEM class
    
    @Attribute len: the length of the LMEM
    @Attribute occur: A list of all occurances in the reference genome, 
                      recorded by the start index in the genome.
    """
    def __init__(self):
        self.len = 0
        self.occur = []

def BWT_search(R, G, start, stop) -> LMEM:
    """
    BWT_search algorithm
    
    @param R: one read
    @param G: reference genome
    @param start: The starting index of the read to perform BWT_search
    @param stop: The stopping index of the read to perform BWT_search
    
    Return: All LMEMs
    """
    pass

def LMEM_explore(R, G, k):
    """
    Find all simple pairs.
    
    @param R: one read
    @param G: reference genome
    @param k: threshold length of a LMEM
    
    Return: list of all LMEMs
    """
    start = 0
    stop = len(R)
    A = []
    for _ in range(stop):
        if R[start] not in "ACGT":
            start += 1
        else:
            LMEM = BWT_search(G, R, start, stop)
            if LMEM.len >= k:
                for p in LMEM.occur:
                    A.append((start, start+LMEM.len-1, p, p+LMEM.len-1, LMEM.len))
                start += LMEM.len
            else:
                start += 1
    return A

def form_cluster(A, g):
    """
    Iterate through all the simple pairs and cluster them based on delta_pos. This finds all candidation
    region alignments and prepare for the normal pair identification.
    
    @param A: List of all LMEMs
    @param g: Threshold for delta_pos
    
    Return: A list of lists. Each of the sublist is a cluster of simple pairs that are g-close to each other.
    """
    A.sort(key=lambda x: x[2]-x[0])
    clustered_A = [[A[0]]]
    
    for i in range(len(A)-1):
        this_sp, next_sp = A[i], A[i+1]
        if (next_sp[2] - next_sp[0]) - (this_sp[2] - this_sp[0]) <= g:
            clustered_A[-1].append(next_sp)
        else:
            clustered_A.append([next_sp])
            
    return clustered_A

def form_candidate_region(clustered_A):
    """
    Remove overlaps between simple pairs and insert normal pairs to form a candidate region.
    
    @param clustered_A: A list of lists. Each of the sublist is a cluster of simple pairs that are 
                        g-close to each other.
                        
    Return: A list of lists. Each sublist is a candidate region.
    """
    
    for cl in clusted_A:
        for i in range(len(cl)-1):
            this_sp, next_sp = cl[i], cl[i+1]
            
            # remove overlap
            if this_sp[3] > next_sp[2]:
                overlap = this_sp[3] - next_sp[2]
                if this_sp[4] < next_sp[4]:
                    this_sp[4] -= overlap
                    this_sp[1] -= overlap
                else:
                    next_sp[3] += overlap
                    this_sp[0] += overlap
            
            # insert normal pair
            if next_sp[0] - this_sp[1] > 1 or next_sp[2] - this_sp[3] > 1:
                normal_pair = (None, None, None, None)
                if next_sp[0] - this_sp[1] > 1:
                    normal_pair[0] = this_sp[1] + 1
                    normal_pair[1] = next_sp[0] - 1
                else:
                    normal_pair[0], normal_pair[1] = -1, -1
                if next_sp[2] - this_sp[3] > 1:
                    normal_pair[2] = this_sp[3] + 1
                    normal_pair[3] = next_sp[2] - 1
                else:
                    normal_pair[2], normal_pair[3] = -1, -1
                    
                cl.insert(i+1,normal_pair)
                
    return clustered_A


            
            
