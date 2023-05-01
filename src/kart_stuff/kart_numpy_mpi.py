from bwa import bwa_run, get_n_reads, get_ref
from copy import deepcopy
from itertools import pairwise

import numpy as np
import random

from mpi4py import MPI

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

def BWT_search(R, G, start, stop, k, verbose=0) -> LMEM:
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
        if verbose == 2:
            print(">> Finding LMEM (Remaining String: {:.2%})".format(end/stop))
            
        if end-start < k:
            return LMEM(0, [])
        
        occur = bwa_run(R[start:end], G, tolerance=0)
        # print("Start", start, "End", end, "Occur", occur)
        if len(occur) != 0:
            break
        else:
            end -= 1
    return LMEM(end-start, occur)

def LMEM_explore(R, G, k, verbose=0):
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
        if verbose >= 1:
            print("> Finding simple pairs ({:.2%})".format(start/stop))
        
        if abs(start - stop) < k:
            break
        
        if R[start] not in "ACGT":
            start += 1
        else:
            LMEM = BWT_search(R, G, start, stop, k, verbose)
            # print(LMEM.occur)
            if LMEM.len >= k:
                for p in LMEM.occur:
                    A.append(np.array([start, start+LMEM.len-1, p, p+LMEM.len-1, LMEM.len]))
                start += LMEM.len
            else:
                start += 1
    if len(A) == 0:
        print("No simple pair found.")
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
    clustered_A = [[A[0]]]
    
    for i in range(len(A)-1):
        this_sp, next_sp = A[i], A[i+1]
        if abs((next_sp[2] - next_sp[0]) - (this_sp[2] - this_sp[0])) <= g:
            clustered_A[-1].append(next_sp)
        else:
            clustered_A.append([next_sp])
    return clustered_A

def form_candidate_regions(clustered_A):
    """
    Remove overlaps between simple pairs and insert normal pairs to form a candidate region.
    
    @param clustered_A: A list of lists. Each of the sublist is a cluster of simple pairs that are 
                        g-close to each other.
                        
    Return: A list of lists. Each sublist is a candidate region.
    """
    
    for idx, cl in enumerate(clustered_A):
        cl.sort(key=lambda x: x[0])
        new_cl = deepcopy(cl)
        
        cnt = 0
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
                normal_pair = np.zeros(4)
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
                
                new_cl.insert(i+1+cnt,normal_pair)
                cnt += 1
                
        clustered_A[idx] = new_cl
    
    return clustered_A

def select_cr(crs, read, ref):
    """
    Select the most probable candidate region from all the candidates
    
    @param crs: List of candidate regions
    @param read: The read that is to be mapped
    @param ref: Reference genome
    
    Return: Two indices, which marks the most probable candidation region.
    """
    best_score = np.NINF
    for cr in crs:
        score = 0
        for p in cr:
            if len(p) == 5:
                score += p[4]
            else:
                score += max(p[0]-p[1], p[2]-p[3]) # penalty for mismatch
        if best_score <= score:
            best_score = score
            best_cr = cr
    
    return best_cr[0][2]-best_cr[0][0], best_cr[-1][3]+len(read)-1-best_cr[-1][1]

def assemble_genome(mapped_list, reads, ref):
    """
    Assemble all the mapped reads.
    
    @param mapped_list: List of all the mapped reads in the form of two indices and the read index
    @param ref: Reference genome
    
    Return: A string. The assembled reads.
    """
    genome = []
    if mapped_list[0][1] > 0:
        genome.append("-"*mapped_list[0][1])
    for this_read, next_read in pairwise(mapped_list):
        read_str = reads[this_read[0]]
        if this_read[2] >= next_read[1]:
            read_str = read_str[:next_read[1]-this_read[1]]
            genome.append(read_str)
        elif next_read[1] - this_read[2] > 1:
            genome.append(read_str)
            genome.append("-"*(next_read[1] - this_read[2] - 1))
    if mapped_list[-1][2] < len(ref)-1:
        genome.append("-"*(len(ref)-1-mapped_list[-1][2]-1))
    
    return "".join(genome)

def kart(read_path, ref_path, k, g, num_reads=100, verbose=0):
    """
    The main function to run the Kart algorithm.
    
    @param read_path: The path to the read file.
    @param ref_path: The path to the reference genome file.
    @param k: Threshold length of a LMEM
    @param g: Threshold for delta_pos
    @param num_read: Number of reads to load for this iteration
    
    Return: An assembled genome.
    
    """
    reads = get_n_reads(read_path, num_reads=num_reads)
    ref = get_ref(ref_path)
    
    mapped_list = []
    
    for read_idx, read in enumerate(reads):
        A = LMEM_explore(read, ref, k, verbose)
        if len(A) == 0:
            print(f"WARNING: Didn't find any simple pair for read {read_idx}")
            print(f"WARNING: Please check read:", reads[read_idx])
            continue
        clustered_A = form_cluster(A, g)
        crs = form_candidate_regions(clustered_A)
        cr = select_cr(crs, read, ref)
        mapped_list.append((read_idx, *cr))
    
    mapped_list.sort(key=lambda x: x[1]) 
    
    assembled = assemble_genome(mapped_list, reads, ref) 
    
    # return assembled

def kart_exp(all_reads, ref_path, k, g, num_reads, verbose=0):
    """
    The main function to run the Kart algorithm.
    
    @param reads: List of reads.
    @param ref_path: The path to the reference genome file.
    @param k: Threshold length of a LMEM
    @param g: Threshold for delta_pos
    @param num_read: Number of reads to load for this iteration
    
    Return: An assembled genome.
    
    """
    ref = get_ref(ref_path)
    
    mapped_list = []
    reads = []
    
    for _ in range(num_reads):
        read = random.choice(all_reads)
        reads.append(read)
    
    for read_idx, read in enumerate(reads):
        A = LMEM_explore(read, ref, k, verbose)
        if len(A) == 0:
            print(f"WARNING: Didn't find any simple pair for read {read_idx}")
            print(f"WARNING: Please check read:", reads[read_idx])
            continue
        clustered_A = form_cluster(A, g)
        crs = form_candidate_regions(clustered_A)
        if verbose != 0:
            pprint_align(crs, read, ref)
        cr = select_cr(crs, read, ref)
        mapped_list.append((read_idx, *cr))
    
    mapped_list.sort(key=lambda x: x[1]) # TODO sort the list based on the reference genome position
    
    assembled = assemble_genome(mapped_list, reads, ref) # TODO assemble_genome function
    
    return assembled

def split_reads(reads, num_workers):
    """
    Split reads into subchunks.
    
    @param reads: a list of all the reads to be mapped
    @param num_workers: number of splits for mpi
    """
    reads_splits = []
    for start, stop in pairwise(np.linspace(0, len(reads), num_workers+1, dtype=int)):
        reads_splits.append(reads[start:stop])
    return reads_splits

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
                else:
                    G_str += G[p[2]:p[3]+1]
                    R_str += R[p[0]:p[1]+1]
                    aligned_str += " " * (len(G[p[2]:p[3]+1]))
        
        G_str += G[cl[-1][3]+1: cl[-1][3]+len(R[cl[-1][1]:])]
        R_str += R[cl[-1][1]+1:]
        
        print("Index:    ", cl[0][2] - cl[0][0], " "*(len(G_str)-4), cl[-1][3]+1+len(R[cl[-1][1]:]))
        print("Reference:", G_str)
        print("Read:     ", R_str)
        print("Aligned:  ", aligned_str)
        
        print()


