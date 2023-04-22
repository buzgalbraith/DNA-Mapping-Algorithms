from Bio import SeqIO
import numpy as np
from numpy import array
import screed
from sourmash import MinHash
from sourmash.minhash import hash_murmur

def split_fastq(fastq,processes):
    fq_dict = SeqIO.index(fastq, "fastq")
    num_keys = len(fq_dict)
    keys_indices_to_keep = getrand_sample(np.arange(num_keys),-(num_keys//-3))#split reads by 1/3rd cus i think its a safe numb
    fq_keys_to_keep = array(list(fq_dict.keys()))[keys_indices_to_keep]
    return fq_keys_to_keep,fq_dict

def build_kmer_index(sequence, ksize):
    n_kmers = len(sequence) - ksize + 1
    mh = MinHash(n=0, ksize=ksize, scaled=1)
    kmer_index = {}
    for i in range(n_kmers):
        if kmer_index.get(hash_murmur(sequence[i:i + ksize])) == None:
            kmer_index[hash_murmur(sequence[i:i + ksize])] = [i]
        else:
            kmer_index[hash_murmur(sequence[i:i + ksize])].append(i)
    return kmer_index,mh

def getrand_sample(indices_array,size):
    return np.random.choice(indices_array,size)

def read_kmers_from_file(filename, ksize):
    kmerized_contigs = [] #each item holds a dictionary of indexed kmers for each contig
    for record in screed.open(filename):
        sequence = record.sequence

        [kmer_index,mh] = build_kmer_index(sequence, ksize)
        kmerized_contigs += [[kmer_index,mh]]

    return kmerized_contigs

def seed_frm_fastq(kmerized_contigs,fq_keys_to_keep,fq_dict,rev_comp,stop_extend):
    read_loc = {}
    contig = kmerized_contigs[0][0]
    mh = kmerized_contigs[0][1]

    for key in fq_keys_to_keep:
        seq_read = fq_dict[key].seq
        if rev_comp:
            seq_read = seq_read.reverse_complement()
        print(key)
        hashed_kmerized_seq_read = mh.seq_to_hashes(str(seq_read),force=True)
        if len(hashed_kmerized_seq_read) < 1: continue
        if contig.get(hashed_kmerized_seq_read[0]) != None:
            i=0
            initial_read_indices = contig[hashed_kmerized_seq_read[0]]

            while i < (len(hashed_kmerized_seq_read)-1) and i <= stop_extend: #by putting this while loop in the extend function, i can maybe minimize repetitive function calls
                print(i,(len(hashed_kmerized_seq_read)-1))
                if contig.get(hashed_kmerized_seq_read[i]) == None:
                    break
                else:
                    place = extend(initial_read_indices,contig,hashed_kmerized_seq_read,i)
                    if place[0]:
                        i+=1
                        if place[1] == 0: loc = place[2]
                    else:
                        place == False
                        break
            if place[0]:
                read_loc[seq_read] = loc
    return read_loc

def extend(initial_read_indices,contig,hashed_kmerized_seq_read,i):
    print('prob line 66')
    if i != 0:
        initial_read_indices = contig[hashed_kmerized_seq_read[i]]
    print('prob line 69')
    if contig.get(hashed_kmerized_seq_read[i+1]) == None:
        place = False
        return [place]
    read_indices_slide = contig[hashed_kmerized_seq_read[i+1]]
    print('prev line fine')
    read_indices_slide = np.asarray(read_indices_slide)
    for loc in initial_read_indices: #the larger the k-mer_size specified the less this loop has to iterate
        idx = (np.abs(read_indices_slide - loc)).argmin()
        if read_indices_slide[idx] - loc == 1:
            place = True
            print('working4')
            return [place,i,loc]
        else:
            place = False
    return [place] #return is outside else statement to allow for loop to find proper index

def switch_dict_keys_values(read_loc): #this is hear b/c of my lack of foresight
    #this switches the items in read_loc, so the location for a read is now the key
    new_read_loc = {}
    for k, v in read_loc.items():
        if v in new_read_loc.keys():
            new_read_loc[v].append(k)
        else:
            new_read_loc[v] = [k]
    return new_read_loc

def create_contig_from_read_loc(new_read_loc):
    sorted_keys = list(new_read_loc.keys())
    sorted_keys.sort()
    reads_to_combine = [str(max(new_read_loc[key])) for key in sorted_keys]
    contig = ''.join(reads_to_combine)
    return contig



