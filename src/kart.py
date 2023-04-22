from bwa import bwa

class LMEM:
    """
    An LMEM class
    
    @Attribute len: the length of the LMEM
    @Attribute occur: A list of all occurances in the reference genome, 
                      recorded by the start index in the genome.
    """
    def __init__(self, length, occur):
        self.len = length
        self.occur = occur

def BWT_search(R, G, start, stop) -> LMEM:
    """
    BWT_search algorithm
    
    @param R: one read
    @param G: reference genome
    @param start: The starting index of the read to perform BWT_search
    @param stop: The stopping index of the read to perform BWT_search
    
    Return: All LMEMs
    """
    end = stop
    while True:
        occur = bwa(R[start:end], G, tolerance=0)
        # print("Start", start, "End", end, "Occur", occur)
        if len(occur) != 0:
            break
        else:
            end -= 1
    return LMEM(end-start, occur)

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
        if start == stop:
            break
        
        if R[start] not in "ACGT":
            start += 1
        else:
            LMEM = BWT_search(R, G, start, stop)
            # print(LMEM.occur)
            if LMEM.len >= k:
                for p in LMEM.occur:
                    A.append([start, start+LMEM.len-1, p, p+LMEM.len-1, LMEM.len])
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
    A.sort(key=lambda x: x[2]-x[0], reverse=True)
    print(A)
    clustered_A = [[A[0]]]
    
    for i in range(len(A)-1):
        this_sp, next_sp = A[i], A[i+1]
        if abs((next_sp[2] - next_sp[0]) - (this_sp[2] - this_sp[0])) <= g:
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
    
    for cl in clustered_A:
        for i in range(len(cl)-1):
            this_sp, next_sp = cl[i], cl[i+1]
            
            # remove overlap
            if this_sp[3] > next_sp[2]:
                overlap = this_sp[3] - next_sp[2]
                if this_sp[4] < next_sp[4]:
                    this_sp[4] -= overlap
                    this_sp[1] -= overlap
                    this_sp[3] -= overlap
                else:
                    next_sp[4] -= overlap
                    next_sp[2] += overlap
                    next_sp[0] += overlap
            
            # insert normal pair
            if next_sp[0] - this_sp[1] > 1 or next_sp[2] - this_sp[3] > 1:
                normal_pair = [None, None, None, None]
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


def pprint_align(clustered_A, R, G):
    """
    Pretty print of candidate alignments
    
    @param clustered_A: candidate regions
    @param R: a read
    @param G: reference genome
    
    Return: None
    """
    for i, cl in enumerate(clustered_A):
        
        print(f"Candidate {i+1}:")
        
        # start = cl[0][2] - cl[0][0]
        # end = cl[-1][3] + (len(R) - cl[-1][1] + 1)
        # print(G[start : end])
        
        G_str = ""
        R_str = ""
        aligned_str = ""
        
        if cl[0][0] != 0:
            G_str += G[cl[0][2] - cl[0][0]: cl[0][2]]
            R_str += R[:cl[0][0]]
            aligned_str += " " * cl[0][0]
            
        for p in cl:
            if len(p) == 5:
                G_str += G[p[2]:p[3]+1]
                R_str += R[p[0]:p[1]+1]
                aligned_str += "*" * (p[4])
            else:
                if p[0] == -1 and p[1] == -1:
                    G_str += G[p[2]:p[3]+1]
                    R_str += "-" * (p[3]-p[2]+1)
                    aligned_str += " " * (p[3]-p[2]+1)
                elif p[2] == -1 and p[3] == -1:
                    G_str += "-"*(p[1]-p[0]+1)
                    R_str += R[p[0]:p[1]+1]
                    aligned_str += " " * (p[1]-p[0]+1)
        
        G_str += G[cl[-1][3]+1: cl[-1][3]+1+len(R[cl[-1][1]:])]
        R_str += R[cl[-1][1]:]
        
        print("Index:    ", cl[0][2] - cl[0][0], " "*(len(G_str)-4), cl[-1][3]+1+len(R[cl[-1][1]:]))
        print("Reference:", G_str)
        print("Read:     ", R_str)
        print("Aligned:  ", aligned_str)
        
        print()
        
            
