from bwa import bwa_run, get_n_reads, get_ref
from copy import deepcopy
from itertools import pairwise

from functools import partial
from multiprocessing.pool import Pool

import numpy as np
import random

from kart_numpy import LMEM_explore, form_cluster, form_candidate_regions, pprint_align, select_cr, assemble_genome

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
    
    LMEM_explore_mp = partial(LMEM_explore, G=ref, k=k, verbose=verbose)
    
    for _ in range(num_reads):
        read = random.choice(all_reads)
        reads.append(read)
    
    with Pool(6) as p:
        As = p.map(LMEM_explore_mp, reads)

    for read_idx, A in enumerate(As):
        if len(A) == 0:
            print(f"WARNING: Didn't find any simple pair.")
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
        
            
