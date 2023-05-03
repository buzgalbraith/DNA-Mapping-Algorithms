import numpy as np
from mpi4py import MPI
import kmer_readmapping as k
import collections
from time import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def spltcontigs(filename,ksize,processes):
    kmerized_contigs = k.read_kmers_from_file(filename,ksize)
    hash_keys,n = np.array(list(kmerized_contigs[0][0].keys()),dtype='uint64').reshape(processes,int(len(kmerized_contigs[0][0].keys())/processes)),len(kmerized_contigs[0][0].keys())
    return [kmerized_contigs,hash_keys,n]

filename,ksize,p = 'SRR11528307_SarsCov2.fasta',15,size #working with 5 processes, cus the contigs kmers divide nicely
contig_mh_hash_keys_n = spltcontigs(filename,ksize,p)
n = contig_mh_hash_keys_n[2]

if rank == 0:
    out = spltcontigs(filename,ksize,p)
    contig_mh, hash_keys, n = out[0],out[1],out[2]
else:
    hash_keys = None

recvbuf = np.empty(int(n/p),dtype='uint64')

print(f'Sending to {p} processes')
comm.Scatter(hash_keys,recvbuf, root = 0)

def call_mappingfunctions(recvbuf,contig_mh_hash_keys_n):
    fastq = 'SRR11528307_R1.fastq'
    fq_dict = k.load_fastq(fastq)
    stop_extend ,rev_comp = 5,False
    kmerized_contig = contig_mh_hash_keys_n[0]
    contig = {key:kmerized_contig[0][0][key] for key in recvbuf}
    kmerized_contig_p = [[contig,kmerized_contig[0][1]],None]
    print('Seeding and extending..')
    read_loc = k.seed_frm_fastq(kmerized_contig_p,fq_dict,rev_comp)
    new_read_loc = k.switch_dict_keys_values(read_loc)
    return new_read_loc

ts = time()
new_read_loc = call_mappingfunctions(recvbuf,contig_mh_hash_keys_n)
new_read_loc_l = np.array(new_read_loc)

print('Gathering data from processes..')
read_locations = comm.gather(new_read_loc,root=0)

if rank ==0:
    collects_dicts = collections.defaultdict(list)
    for dict in read_locations:
        for key,v in dict.items():
            collects_dicts[key].append(v)
    print('Creating contig..')
    contig = k.fltr_found_reads_makeContig(collects_dicts, filename)

    print(' Took {}s'.format(time() - ts))

quit()


