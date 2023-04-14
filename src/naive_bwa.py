def bwa_run(band_start, band_end, read, index, bwt_reference, bands):
    """
    runs the bwa algorithm on a single read. 
    args"
    band_start, band_end: (int) index postions of band in reference genome 
    read: (str)  a single sequence produced from a sequencer.
    index: (int) current position in read 
    """

    if index == -1:
        return band_start, band_end 
    if band_start == band_end:
        return -1
    s = read[index]
    rank_top = rank(bwt_ref,s,band_start)
    rank_bottom = rank(bwt_ref,s,band_end)
    return bwa_run(bands[read[index]]+rank_top, bands[read[index]]+rank_bottom, read , index-1, bands=bands,bwt_reference=bwt_reference)


def build_bands(ref):
    """
    builds the bands for a given reference genome 
    args:
    ref: (str) refereence genome
    """
    sorted_ref = sorted(ref)
    bands = {}
    prev_char = None
    prev_index = 0
    for i, char in enumerate(sorted_ref):
        if char != prev_char:
            bands[prev_char] = prev_index+1
            prev_char = char
            prev_index = i
    bands[prev_char] = prev_index
    return bands


def suffix_array(reference, v='$'):
    """
    builds the suffix array 
    args:
    ref: (str) refereence genome
    v: (char) charicter not in refernce array to be added. defaults to $
    """
    
    assert v not in reference,  "{0} already in input string".format(v)
    reference = reference+v
    suffixes = [reference[i:] for i in range(len(reference))]
    suffix_array = sorted(range(len(reference)), key=lambda i: suffixes[i])
    bwt = ''.join([reference[i-1] for i in suffix_array])
    return bwt, suffix_array

def bwa(read, reference, v="$"):
     """
    whole bwa algorithm
    args:
    read: (str) indivudal reading
    ref: (str) refereence genome
    v: (char) charicter not in refernce array to be added. defaults to $
    """
    bands = build_band_array(read,reference)
    backwards_read = read[::-1]
    bwt_ref, suffix_ref=suffix_array(reference,v=v)
    indecies= bwa_run(bands[backwards_read[0]], bands[backwards_read[1]], read, index=len(read)-2, bwt_reference=bwt_ref , bands=bands)
    if indecies==-1:
        return "No match"
    return suffix_ref[indecies[0]]