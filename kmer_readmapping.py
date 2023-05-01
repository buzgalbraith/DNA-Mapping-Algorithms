from Bio import SeqIO
import numpy as np
import screed
from sourmash import MinHash
from sourmash.minhash import hash_murmur

def load_fastq(fastq):
    fq_dict = SeqIO.index(fastq, "fastq")
    return fq_dict

def build_kmer_index(sequence, ksize): #index, create kmers, and hash kmers from the sequence
    n_kmers = len(sequence) - ksize + 1
    mh = MinHash(n=0, ksize=ksize, scaled=1)
    kmer_index = {}
    for i in range(n_kmers):
        if kmer_index.get(hash_murmur(sequence[i:i + ksize])) == None:
            kmer_index[hash_murmur(sequence[i:i + ksize])] = [i]
        else:
            kmer_index[hash_murmur(sequence[i:i + ksize])].append(i)
    return kmer_index,mh

def read_kmers_from_file(filename, ksize):
    kmerized_contigs = [] #each item holds a dictionary of indexed kmers for each contig
    for record in screed.open(filename):
        sequence = record.sequence

        [kmer_index,mh] = build_kmer_index(sequence, ksize)
        kmerized_contigs += [[kmer_index,mh]]

    return kmerized_contigs

def seed_frm_fastq(kmerized_contigs,fq_dict,rev_comp):
    read_loc = {}
    contig = kmerized_contigs[0][0] #only working with first contig for this implementation
    mh = kmerized_contigs[0][1]

    for key in fq_dict.keys():
        seq_read = fq_dict[key].seq
        if rev_comp:
            seq_read = seq_read.reverse_complement()
        #print(key)
        hashed_kmerized_seq_read = mh.seq_to_hashes(str(seq_read),force=True)
        if len(hashed_kmerized_seq_read) < 1: continue
        if contig.get(hashed_kmerized_seq_read[0]) != None:
            i=0
            initial_read_indices = contig[hashed_kmerized_seq_read[0]]

            while i < (len(hashed_kmerized_seq_read)-1): #by putting this while loop in the extend function, i can maybe minimize repetitive function calls
                #print(i,(len(hashed_kmerized_seq_read)-1))
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
    if i != 0:
        initial_read_indices = contig[hashed_kmerized_seq_read[i]]
    if contig.get(hashed_kmerized_seq_read[i+1]) == None:
        place = False
        return [place]

    read_indices_slide = contig[hashed_kmerized_seq_read[i+1]]
    read_indices_slide = np.asarray(read_indices_slide)
    for loc in initial_read_indices: #the larger the k-mer_size specified the less this loop has to iterate
        idx = (np.abs(read_indices_slide - loc)).argmin()
        if read_indices_slide[idx] - loc == 1:
            place = True
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

def fltr_found_reads_makeContig(new_read_loc,filename):
    sorted_keys = list(new_read_loc.keys())
    sorted_keys.sort()
    reads_to_combine = []
    i = 0
    while i < len(sorted_keys):
        key = sorted_keys[i]
        read = max(new_read_loc[key]) #grab largest read for a particular index
        if i > 0 and key >= len(reads_to_combine[-1]) + tmp: #make sure reads dont overlap when combining
            num_gaps = key - tmp -1 - len(reads_to_combine[-1])
            reads_to_combine.extend(['-'] * num_gaps)
            reads_to_combine.append(str(read))
            tmp = key
        if i == 0:
            num_gaps = key - 1 #gaps between first nucleotide and where first contig is found
            reads_to_combine.extend(['-']*num_gaps)
            reads_to_combine.append(str(read))
            tmp = key
        i += 1

    fasta = [record.sequence for record in screed.open(filename)]
    true_contig_size = len(fasta[0])

    out_contig = ''.join(reads_to_combine)
    gaps_tofill = true_contig_size - len(out_contig)
    gaps = '-' * gaps_tofill
    res = out_contig + gaps
    return res


def run(): #runs functions in proper order
    fasta,fastq = 'SRR11528307_SarsCov2.fasta','SRR11528307_R1.fastq'
    ksize = 15
    rev_comp = False

    fq_dict = load_fastq(fastq)
    kmerized_contigs = read_kmers_from_file(fasta, ksize)
    print('Indexed contig...')
    read_loc = seed_frm_fastq(kmerized_contigs,fq_dict,rev_comp)
    print('Read kmers seeded and extended..')
    read_loc_corrected = switch_dict_keys_values(read_loc)
    contig = fltr_found_reads_makeContig(read_loc_corrected,fasta)
    print('Combined found reads into a contig')

    return contig

#contig = run()


