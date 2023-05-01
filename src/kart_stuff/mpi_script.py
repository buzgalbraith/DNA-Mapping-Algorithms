
from kart_numpy_mpi import split_reads, LMEM_explore, form_cluster, form_candidate_regions, pprint_align, select_cr, assemble_genome
from bwa import get_ref
import random
import pickle
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
num_workers = comm.Get_size()
rank = comm.Get_rank()

if rank==0:
    with open("100_reads_perturbed.pt", "rb") as f:    
        reads = pickle.load(f)
    start = time.time()
    read_data = split_reads(reads[:20], num_workers) # TODO split_reads into sublists

else:
    read_data = None 

read_data_chunk = comm.scatter(read_data, root=0)

ref = get_ref(r'data/SRR11528307_SarsCov2.fna')
for read in read_data_chunk:
    A = LMEM_explore(read, ref, k=10, verbose=0)
    clustered_A = form_cluster(A, g=10)
    crs = form_candidate_regions(clustered_A)
    cr = select_cr(crs, read, ref)
mapped_list = comm.gather(cr, root=0)

if rank==0:
    print(mapped_list)
    end = time.time()
    print("Time spent:", end-start)
